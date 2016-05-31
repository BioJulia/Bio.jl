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
An alignment data for the BAM file format.
"""
type BAMAlignment <: IntervalTrees.AbstractInterval{Int64}
    seqname::StringField
    position::Int64
    mapq::UInt8
    bin::UInt16
    flag::UInt16
    next_seqname::StringField
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
end

function BAMAlignment()
    return BAMAlignment(
        StringField(), -1, 0, 0, 0,
        StringField(), -1, -1,
        UInt8[], -1, -1, -1, -1, -1)
end

function Base.show(io::IO, aln::BAMAlignment)
    println(io, summary(aln), ':')
    println(io, "  seqname:         ", aln.seqname)
    println(io, "  position:        ", aln.position)
    println(io, "  mapping quality: ", aln.mapq)
    println(io, "  bin:             ", aln.bin)
    println(io, "  flag:            ", bin(aln.flag, 16))
    println(io, "  next seqname:    ", aln.next_seqname)
    println(io, "  next position:   ", aln.next_pos)
    println(io, "  template length: ", aln.tlen)
    println(io, "  sequence:        ", string(sequence(aln)))
    println(io, "  CIGAR:           ", cigar(aln))
    println(io, "  quality scores:  ", qualities(aln))
      print(io, "  auxiliary data:  ", auxiliary(aln))
end

# Interval interface
function Base.first(bam::BAMAlignment)
    return bam.position
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

"""
Two-byte tag of auxiliary data in SAM and BAM.
"""
immutable AuxTag
    data::Tuple{UInt8,UInt8}
end

AuxTag(x::UInt8, y::UInt8) = AuxTag((x, y))
AuxTag(x::Char, y::Char) = AuxTag(UInt8(x), UInt8(y))

function Base.getindex(tag::AuxTag, i::Integer)
    if i == 1
        return tag.data[1]
    elseif i == 2
        return tag.data[2]
    end
    throw(BoundsError(i))
end

function Base.show(io::IO, tag::AuxTag)
    write(io, '"', tag[1], tag[2], '"')
    return
end

immutable AuxDataDict <: Associative{AuxTag,Any}
    data::Vector{UInt8}
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

function Base.getindex(dict::AuxDataDict, tag::AbstractString)
    if length(tag) != 2
        error("BAM auxillary data tags must be of length 2")
    end
    return _auxiliary(dict.data, 1, UInt8(tag[1]), UInt8(tag[2]))
end

function Base.getindex(dict::AuxDataDict, key::AuxTag)
    return auxiliary(dict.data, key[1], key[2])
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

function auxiliary(aln::BAMAlignment, tag::AbstractString)
    if length(tag) != 2
        error("BAM auxillary data tags must be of length 2")
    end
    return auxiliary(aln, UInt8(tag[1]), UInt8(tag[2]))
end

# Return a specific tag
function auxiliary(aln::BAMAlignment, t1::UInt8, t2::UInt8)
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

immutable BAMParser{T<:BufferedInputStream} <: AbstractParser
    stream::T
    header_text::StringField
    refs::Vector{Tuple{StringField,Int}}
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
    text = StringField(read!(stream, Vector{UInt8}(textlen)))

    # reference sequences
    n_refs = read(stream, Int32)
    refs = Vector{Tuple{StringField,Int}}()
    for i in 1:n_refs
        namelen = read(stream, Int32)
        name = StringField(read!(stream, Vector{UInt8}(namelen)))
        reflen = read(stream, Int32)
        push!(refs, (name, reflen))
    end

    return BAMParser(stream, text, refs)
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

    aln.position       = fields.pos + 1  # make 1-based
    aln.mapq           = (fields.bin_mq_nl >> 8) & 0xff
    aln.bin            = fields.bin_mq_nl >> 16
    aln.flag           = fields.flag_nc >> 16
    aln.next_pos       = fields.next_pos
    aln.tlen           = fields.tlen
    aln.cigar_position = l_read_name + 1  # read name is NULL-terminated
    aln.seq_position   = aln.cigar_position + 4n_cigar_op
    aln.seq_length     = fields.l_seq
    aln.qual_position  = aln.seq_position + cld(fields.l_seq, 2)
    aln.aux_position   = aln.qual_position + fields.l_seq

    if fields.refid ≥ 0
        aln.seqname = parser.refs[fields.refid+1][1]
    else
        empty!(aln.seqname)
    end
    if fields.next_refid ≥ 0
        aln.next_seqname = parser.refs[fields.next_refid+1][1]
    else
        empty!(aln.next_seqname)
    end

    return aln
end


immutable BAMWriter{T<:BufferedInputStream}
    stream::T
    header_text::StringField
    refs::Vector{Tuple{StringField,Int}}
end

# write a BAM header to `writer.stream`.
function write_header(writer::BAMWriter)
    stream = writer.stream

    # magic
    write(stream, "BAM\0")

    # header
    write(stream, Int32(length(writer.header_text)))
    write(stream, writer.header_text)

    # reference sequences
    write(stream, Int32(length(writer.refs)))
    for (refname, reflen) in writer.refs
        write(stream, Int32(length(refname)))
        write(stream, refname)
        write(stream, Int32(reflen))
    end
end

function Base.write(writer::BAMWriter, aln::BAMAlignment)
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
    n += write(stream, Int32(aln.position))
    n += write(stream, UInt32(bin_mq_nl))
    n += write(stream, UInt32(flag_nc))
    n += write(stream, Int32(aln.seq_length))
    n += write(stream, Int32(aln.next_refid))
    n += write(stream, Int32(aln.next_pos))
    n += write(stream, Int32(aln.tlen))
    n += write(stream, aln.data)

    return n
end
