# BAM
# ===
#
# The BAM file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
The BAM file format.

See "Sequence Alignment/Map Format Specification" for details:
<http://samtools.github.io/hts-specs/SAMv1.pdf>
"""
immutable BAM <: FileFormat end


"""
Two-byte tag of data used in SAM and BAM.
"""
immutable Tag <: AbstractString
    data::Tuple{UInt8,UInt8}
end

function Tag(x::UInt8, y::UInt8)
    return Tag((x, y))
end

function Base.convert(::Type{Tag}, s::AbstractString)
    if length(s) != 2
        throw(ArgumentError("tag must be of length 2"))
    end
    return Tag(UInt8(s[1]), UInt8(s[2]))
end

Base.length(tag::Tag) = 2
Base.endof(tag::Tag) = 2
Base.start(tag::Tag) = 1
Base.done(tag::Tag, i::Int) = i > 2
Base.next(tag::Tag, i::Int) = tag.data[i], i + 1

function twobytes(tag::Tag)
    return tag.data[1], tag.data[2]
end

function Base.checkbounds(tag::Tag, i::Integer)
    if 1 ≤ i ≤ endof(tag)
        return true
    end
    throw(BoundsError(i))
end

function Base.getindex(tag::Tag, i::Int)
    checkbounds(tag, i)
    return tag.data[i]
end

function Base.write(io::IO, tag::Tag)
    return write(io, tag[1], tag[2])
end


"""
A header for SAM/BAM file format.
"""
type SAMHeader <: Associative{Tag,Any}
    # expected value types for each tag:
    #   * @HD:    Associative{Tag,Any}
    #   * @CO:    Vector{ASCIIString}
    #   * others: Vector{Associative{Tag,Any}}
    dict::OrderedDict{Tag,Any}
end

function Base.keys(header::SAMHeader)
    return keys(header.dict)
end

function Base.getindex(header::SAMHeader, tag)
    return getindex(header.dict, convert(Tag, tag))
end

function Base.setindex!(header::SAMHeader, value, tag)
    return setidnex!(header.dict, value, convert(Tag, tag))
end

function Base.length(header::SAMHeader)
    return length(header.dict)
end

# parse SAM header and return the result as a nested dictionary.
function parse_samheader(io::IO)
    # keep the order of header lines
    d = OrderedDict{Tag,Any}()
    while !eof(io)
        mark(io)
        line::ASCIIString = readline(io)
        if line[1] != '@'
            reset(io)
            break
        end
        tag = line[2:3]
        if line[4] != '\t'
            error("invalid SAM header")
        end
        rest = chomp(line[5:end])
        if tag == "HD"
            d[tag] = parse_samheader_values(rest)
        else
            if !haskey(d, tag)
                d[tag] = []
            end
            if tag == "CO"
                push!(d[tag], rest)
            else
                push!(d[tag], parse_samheader_values(rest))
            end
        end
    end
    return SAMHeader(d)
end

# parse a header line after "@<tag>\t"
function parse_samheader_values(line)
    # keep the order of values
    ret = OrderedDict{Tag,ASCIIString}()
    for pair in split(line, '\t')
        tag = pair[1:2]
        if pair[3] != ':'
            error("invalid SAM header")
        end
        ret[tag] = pair[4:end]
    end
    return ret
end

function check_samheader(header::SAMHeader)
    # TODO
end

# write the SAM format header to `io`.
function Base.write(io::IO, header::SAMHeader)
    check_samheader(header)
    n = 0

    # HD is the first line if present
    if haskey(header, "HD")
        n += write(io, "@HD")
        for (key, val) in header["HD"]
            n += write(io, '\t')
            n += write(io, key, ':', val)
        end
        n += write(io, '\n')
    end

    for (tag, records) in header
        if tag == Tag("HD")
            continue
        end
        for record in records
            n += write_samheader_record(io, tag, record)
        end
    end

    return n
end

# write a header record of `tag` to `io`: a single line of a SAM header.
function write_samheader_record(io, tag, record)
    n = 0
    n += write(io, '@', tag)
    for (key, val) in record
        n += write(io, '\t')
        n += write(io, key, ':', val)
    end
    n += write(io, '\n')
    return n
end


"""
A list of reference sequences in a SAM/BAM file.

This supports following two mappings:
* integer index -> (sequence name, sequence length)
* sequence name -> (integer index, sequence length)
"""
immutable ReferenceSequences <: AbstractVector{Tuple{ASCIIString,Int64}}
    names::Vector{ASCIIString}
    name2index::Dict{ASCIIString,Int}
    seqlens::Vector{Int64}
end

function ReferenceSequences()
    return ReferenceSequences([])
end

function ReferenceSequences(iter)
    names = ASCIIString[]
    name2index = Dict{ASCIIString,Int}()
    seqlens = Int64[]
    for (i, (seqname, seqlen)) in enumerate(iter)
        push!(names, seqname)
        push!(seqlens, seqlen)
        name2index[seqname] = i
    end
    return ReferenceSequences(names, name2index, seqlens)
end

function Base.checkbounds(refseqs::ReferenceSequences, index::Integer)
    if 1 ≤ index ≤ endof(refseqs)
        return true
    end
    throw(BoundsError(index))
end

function Base.getindex(refseqs::ReferenceSequences, index::Integer)
    checkbounds(refseqs, index)
    name = refseqs.names[index]
    seqlen = refseqs.seqlens[index]
    return name, seqlen
end

function Base.getindex(refseqs::ReferenceSequences, name::AbstractString)
    index = refseqs.name2index[name]
    seqlen = refseqs.seqlens[index]
    return index, seqlen
end

function Base.push!(refseqs::ReferenceSequences, pair::Tuple{AbstractString,Integer})
    name, seqlen = pair
    push!(refseqs.names, name)
    refseqs.name2index[name] = length(refseqs.names)
    push!(refseqs.seqlens, seqlen)
    return refseqs
end

Base.size(refseqs::ReferenceSequences) = (length(refseqs.names),)


"""
An alignment data for the BAM file format.
"""
type BAMAlignment <: IntervalTrees.AbstractInterval{Int64}
    refid::Int32
    pos::Int64
    mapq::UInt8
    bin::UInt16
    flag::UInt16
    next_refid::Int32
    next_pos::Int64
    tlen::Int32

    # variable length data:
    #   * read name
    #   * cigar data
    #   * sequence
    #   * quality scores
    #   * and aux
    # offsets within the data
    data::Vector{UInt8}
    cigar_position::Int32
    seq_position::Int32
    seq_length::Int32
    qual_position::Int32
    aux_position::Int32

    # reference sequences in BAM file
    refseqs::ReferenceSequences
end

function BAMAlignment()
    return BAMAlignment(
        0, -1, 0, 0, 0,
        0, -1, -1,
        UInt8[], -1, -1, -1, -1, -1,
        ReferenceSequences()
    )
end

function Base.show(io::IO, aln::BAMAlignment)
    if aln.refid > 0
        seqname, = aln.refseqs[aln.refid]
    else
        seqname = ""
    end
    if aln.next_refid > 0
        next_seqname, = aln.refseqs[aln.next_refid]
    else
        next_seqname = ""
    end

    println(io, summary(aln), ':')
    println(io, "  seqname:         ", seqname)
    println(io, "  position:        ", aln.pos)
    println(io, "  mapping quality: ", aln.mapq)
    println(io, "  bin:             ", aln.bin)
    println(io, "  flag:            ", bin(aln.flag, 16))
    println(io, "  next seqname:    ", next_seqname)
    println(io, "  next position:   ", aln.next_pos)
    println(io, "  template length: ", aln.tlen)
    println(io, "  sequence:        ", string(sequence(aln)))
    println(io, "  CIGAR:           ", cigar(aln))
    println(io, "  quality scores:  ", qualities(aln))
      print(io, "  auxiliary data:  ", auxiliary(aln))
end

# Interval interface
function Base.first(aln::BAMAlignment)
    return aln.pos
end

function Base.last(bam::BAMAlignment)
    # TODO: calculate end
end


"""
Return a `StringField` containing the read name. This becomes invalidated if
the `BAMAlignment` is overwritten.
"""
function name(bam::BAMAlignment)
    # '-2' because this is stored null-terminated
    return StringField(bam.data[1:bam.cigar_position-2])
end

"""
Copy the read's name into `name`.
"""
function readname!(bam::BAMAlignment, name::StringField)
    copy!(name, bam.data, 1, bam.cigar_position-1)
end

function cigar(aln::BAMAlignment)
    cigar = CIGAR()
    data = aln.data
    for i in aln.cigar_position:4:aln.seq_position-1
        x = UInt32(data[i  ])       |
            UInt32(data[i+1]) <<  8 |
            UInt32(data[i+2]) << 16 |
            UInt32(data[i+3]) << 24
        len = x >> 4
        op = Operation(x & 0b1111)
        push!(cigar, (op, len))
    end
    return cigar
end

# TODO: DNA_Gap
# "=ACMGRSVTWYHKDBN" -> [0,16)
const bam_nucs = [
    DNA_Gap, DNA_A, DNA_C, DNA_M,
    DNA_G,   DNA_R, DNA_S, DNA_V,
    DNA_T,   DNA_W, DNA_Y, DNA_H,
    DNA_K,   DNA_D, DNA_B, DNA_N
]

function sequence(aln::BAMAlignment)
    return sequence!(aln, DNASequence(aln.seq_length))
end

function sequence!(aln::BAMAlignment, seq::DNASequence)
    if length(seq) != aln.seq_length
        seq = resize!(seq, aln.seq_length)
    end
    data = aln.data
    i = 2
    j = aln.seq_position
    while i ≤ aln.seq_length
        x = data[j]
        seq[i-1] = bam_nucs[(x >>   4) + 1]
        seq[i  ] = bam_nucs[(x & 0x0f) + 1]
        i += 2
        j += 1
    end
    if isodd(aln.seq_length)
        x = data[j]
        seq[i-1] = bam_nucs[(x >>   4) + 1]
    end
    return seq
end

function qualities(aln::BAMAlignment)
    r = aln.qual_position:aln.aux_position-1
    qs = Vector{Int8}(length(r))
    @inbounds for (i, j) in enumerate(r)
        qs[i] = aln.data[j] - 33
    end
    return qs
end

function qualities!(aln::BAMAlignment, qs::Vector{Int8})
    r = aln.qual_position:aln.aux_position-1
    if length(qs) != length(r)
        resize!(qs, length(r))
    end
    @inbounds for (i, j) in enumerate(r)
        qs[i] = aln.data.data[j] - 33
    end
    return qs
end


const auxtype = Dict{UInt8,DataType}(
    'A' => Char,
    'c' => Int8,
    'C' => UInt8,
    's' => Int16,
    'S' => UInt16,
    'i' => Int32,
    'I' => UInt32,
    'f' => Float32,
    'd' => Float64,
    'Z' => ASCIIString
)

"""
Auxiliary data dictionary of SAM/BAM file formats.

This is not designed for very large dictionaries: time complexities are O(N) in lookup and update operations.
"""
immutable AuxDataDict <: Associative{Tag,Any}
    data::Vector{UInt8}
end

function Base.getindex(dict::AuxDataDict, key)
    t1, t2 = twobytes(Tag(key))
    return _auxiliary(dict.data, 1, t1, t2)
end

Base.eltype(::Type{AuxDataDict}) = Tuple{AuxTag,Any}
Base.length(dict::AuxDataDict) = count_auxtags(dict.data, 1)
Base.start(dict::AuxDataDict) = 1
Base.done(dict::AuxDataDict, pos) = pos > length(dict.data)
function Base.next(dict::AuxDataDict, pos)
    data = dict.data
    tag = AuxTag(data[pos], data[pos+1])
    pos, typ = getauxtype(data, pos + 2)
    pos, value = getauxdata(data, pos, typ)
    return (tag, value), pos
end

function auxiliary(aln::BAMAlignment)
    return AuxDataDict(aln.data[aln.aux_position:end])
end

function auxiliary(aln::BAMAlignment, tag)
    t1, t2 = twobytes(Tag(tag))
    return _auxiliary(aln.data, aln.aux_position, t1, t2)
end

function _auxiliary(data::Vector{UInt8}, pos::Integer, t1::UInt8, t2::UInt8)
    p::Int = pos

    while p ≤ length(data) && (data[p] != t1 || data[p+1] != t2)
        p = next_tag_position(data, p)
    end

    if p > length(data)
        throw(KeyError(AuxTag(t1, t2)))
    else
        p, typ = getauxtype(data, p + 2)
        _, value = getauxdata(data, p, typ)
        return value
    end
end

function getauxtype(data::Vector{UInt8}, p::Int)
    t = data[p]
    if t == UInt8('B')
        return p + 2, Vector{auxtype[data[p+1]]}
    else
        return p + 1, auxtype[t]
    end
end

function getauxdata{T}(data::Vector{UInt8}, p::Int, ::Type{T})
    return p + sizeof(T), unsafe_load(Ptr{T}(pointer(data, p)))
end

function getauxdata(data::Vector{UInt8}, p::Int, ::Type{Char})
    return p + 1, Char(unsafe_load(pointer(data, p)))
end

function getauxdata{T}(data::Vector{UInt8}, p::Int, ::Type{Vector{T}})
    n = unsafe_load(Ptr{Int32}(pointer(data, p)))
    p += 4
    xs = Array(T, n)
    unsafe_copy!(pointer(xs), Ptr{T}(pointer(data, p)), n)
    return p + n * sizeof(T), xs
end

function getauxdata(data::Vector{UInt8}, p::Int, ::Type{ASCIIString})
    dataptr = pointer(data, p)
    endptr = ccall(:memchr, Ptr{Void}, (Ptr{Void}, Cint, Csize_t),
                   dataptr, '\0', length(data) - p + 1)
    q = p + (endptr - dataptr) - 1
    return q + 2, ASCIIString(data[p:q])
end

function count_auxtags(data::Vector{UInt8}, p::Int)
    count = 0
    while p ≤ length(data)
        count += 1
        p = next_tag_position(data, p)
    end
    return count
end

# Find the starting position of a next tag in `data` after `p`.
# `(data[p], data[p+1])` is supposed to be a current tag.
function next_tag_position(data::Vector{UInt8}, p::Int)
    typ = Char(data[p+2])
    p += 3
    if typ == 'A'
        p += 1
    elseif typ == 'c' || typ == 'C'
        p += 1
    elseif typ == 's' || typ == 'S'
        p += 2
    elseif typ == 'i' || typ == 'I'
        p += 4
    elseif typ == 'f'
        p += 4
    elseif typ == 'd'
        p += 8
    elseif typ == 'Z' || typ == 'H'
        while data[p] != 0x00  # NULL-terminalted string
            p += 1
        end
        p += 1
    elseif typ == 'B'
        eltyp = Char(data[p])
        elsize = eltyp == 'c' || eltyp == 'C'                 ? 1 :
                 eltyp == 's' || eltyp == 'S'                 ? 2 :
                 eltyp == 'i' || eltye == 'I' || eltyp == 'f' ? 4 :
                 error("unrecognized auxiliary type: ", eltyp)
        p += 1
        n = unsafe_load(Ptr{Int32}(pointer(data, p)))
        p += elsize * n
    else
        error("unrecognized auxiliary type: ", typ)
    end
    return p
end


"""
A parser for BAM file format.
"""
immutable BAMParser{T<:BufferedInputStream} <: AbstractParser
    stream::T
    header::SAMHeader
    refseqs::ReferenceSequences
end

Base.eof(parser::BAMParser) = eof(parser.stream)

function Base.open(source_stream::BufferedInputStream, ::Type{BAM})
    stream = BufferedInputStream(BGZFSource(source_stream))

    # magic bytes
    B = read(stream, UInt8)
    A = read(stream, UInt8)
    M = read(stream, UInt8)
    if B != UInt8('B') || A != UInt8('A') || M != UInt8('M') ||
        read(stream, UInt8) != 0x01
        error("input was not a valid BAM file")
    end

    # header text
    textlen = read(stream, Int32)
    header = parse_samheader(IOBuffer(read(stream, UInt8, textlen)))

    # reference sequences
    refseqs = ReferenceSequences()
    n_refs = read(stream, Int32)
    for _ in 1:n_refs
        namelen = read(stream, Int32)
        # remove the last NULL character
        seqname = chop(bytestring(read(stream, UInt8, namelen)))
        seqlen = read(stream, Int32)
        push!(refseqs, (seqname, seqlen))
    end

    return BAMParser(stream, header, refseqs)
end


# Leading size parts of a bam entry
immutable BAMEntryHead
    refid::Int32
    pos::Int32
    bin_mq_nl::UInt32
    flag_nc::UInt32
    l_seq::Int32
    next_refid::Int32
    next_pos::Int32
    tlen::Int32
end

function Base.read!(parser::BAMParser, aln::BAMAlignment)
    stream = parser.stream
    if eof(stream)
        throw(EOFError())
    end

    block_size = read(stream, Int32)
    if !BufferedStreams.ensurebuffered!(stream, block_size)
        throw(EOFError())
    end

    # this hack significantly improves the performance of reading alignment data
    ptr::Ptr{BAMEntryHead} = pointer(stream.buffer, stream.position)
    fields = unsafe_load(ptr)
    stream.position += sizeof(BAMEntryHead)

    datalen = block_size - sizeof(BAMEntryHead)
    if length(aln.data) < datalen
        resize!(aln.data, datalen)
    end
    n = readbytes!(stream, aln.data, datalen)
    @assert n == datalen

    n_cigar_op = fields.flag_nc & 0xffff
    l_read_name = fields.bin_mq_nl & 0xff

    # TODO: directly map the data buffer to a BAMAlignment may be faster
    aln.refid          = fields.refid + 1  # make 1-based
    aln.pos            = fields.pos + 1  # make 1-based
    aln.mapq           = (fields.bin_mq_nl >> 8) & 0xff
    aln.bin            = fields.bin_mq_nl >> 16
    aln.flag           = fields.flag_nc >> 16
    aln.next_refid     = fields.next_refid + 1
    aln.next_pos       = fields.next_pos + 1
    aln.tlen           = fields.tlen
    aln.cigar_position = l_read_name + 1  # read name is NULL-terminated
    aln.seq_position   = aln.cigar_position + 4n_cigar_op
    aln.seq_length     = fields.l_seq
    aln.qual_position  = aln.seq_position + cld(fields.l_seq, 2)
    aln.aux_position   = aln.qual_position + fields.l_seq
    aln.refseqs        = parser.refseqs

    return aln
end


type BAMWriter{T<:BufferedOutputStream}
    stream::T
    header::SAMHeader
end

function Base.open(output::BufferedOutputStream, ::Type{BAM};
                   header::SAMHeader=SAMHeader())
    stream = BufferedOutputStream(BGZFSink(output))
    writer = BAMWriter(stream, header)
    write_header(stream, header)
    return writer
end

# write a BAM header to `writer.stream`.
function write_header(writer::BAMWriter)
    stream = writer.stream
    n = 0

    # magic
    n += write(stream, "BAM\0")

    # header
    buf = IOBuffer()
    write(buf, writer.header)
    header_array = takebuf_array(buf)
    n += write(stream, Int32(length(header_array)))
    n += write(stream, header_array)

    # reference sequences
    refseqs = writer.header["SQ"]
    n += write(stream, Int32(length(refseqs)))
    for record in refseqs
        name = record["SN"]
        len  = record["LN"]
        n += write(stream, Int32(length(name) + 1))
        # sequence name must be NULL-terminated
        n += write(stream, name, '\0')
        n += write(stream, Int32(len))
    end

    return n
end

function Base.write(writer::BAMWriter, aln::BAMAlignment)
    if !(0 ≤ aln.refid ≤ endof(writer.refseqs))
        error(aln.refid, " is not in the reference sequence set")
    elseif !(0 ≤ aln.next_refid ≤ endof(writer.refseqs))
        error(aln.next_refid, " is not in the reference sequence set")
    end

    stream = writer.stream

    block_size = sizeof(BAMEntryHead) + sizeof(aln.data)
    bin_mq_nl = UInt32(aln.bin)  << 16 |
                UInt32(aln.mapq) <<  8 |
                UInt32(aln.cigar_position - 1)
    @assert rem(aln.seq_position - aln.cigar_position, 4) == 0
    n_cigar_op = div(aln.seq_position - aln.cigar_position, 4)
    flag_nc = UInt32(aln.flag) << 16 | UInt32(n_cigar_op)

    n = 0
    n += write(stream, Int32(block_size))
    n += write(stream, Int32(aln.refid))
    n += write(stream, Int32(aln.pos))
    n += write(stream, UInt32(bin_mq_nl))
    n += write(stream, UInt32(flag_nc))
    n += write(stream, Int32(aln.seq_length))
    n += write(stream, Int32(aln.next_refid))
    n += write(stream, Int32(aln.next_pos))
    n += write(stream, Int32(aln.tlen))
    n += write(stream, aln.data)

    return n
end
