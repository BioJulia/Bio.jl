using Libz
using BufferedStreams

using Bio.Seq
using Bio.Seq: DNA_INVALID
using Bio.Align
#using Bio.BGZF
using Bio.StringFields
using Bio.Intervals
using IntervalTrees
using BGZF

immutable BAM <: Bio.FileFormat end


type BAMAlignment <: AbstractInterval{Int64}
    seqname::StringField
    position::Int64
    mapq::UInt8
    bin::UInt16
    flag::UInt16
    next_seqname::StringField
    next_pos::Int64
    tlen::Int32

    # variable length data: read name, cigar data, sequence, quality scores, and aux
    # offsets within the data
    #data::StringField
    data::Vector{UInt8}
    cigar_position::Int32
    seq_position::Int32
    seq_length::Int32
    qual_position::Int32
    aux_position::Int32
end

function BAMAlignment()
    #return BAMAlignment(StringField(), -1, 0, 0, 0, StringField(), -1, -1,
    #                    StringField(), -1, -1, -1, -1, -1)
    return BAMAlignment(StringField(), -1, 0, 0, 0, StringField(), -1, -1,
                        UInt8[], -1, -1, -1, -1, -1)
end

# Interval interface
function first(bam::BAMAlignment)
    return bam.position
end


function last(bam::BAMAlignment)
    # TODO: calculate end
end


"""
Return a `StringField` containing the read name. This becomes invalidated if
the `BAMAlignment` is overwritten.
"""
#function readname(bam::BAMAlignment)
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

function Align.cigar(aln::BAMAlignment)
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


const bam4bit_to_char = [
    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N' ]


const bam4bit_to_dna = [
    DNA_INVALID, DNA_A,       DNA_C,       DNA_INVALID,
    DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
    DNA_T,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
    DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_N ]

#=
"""
Return a DNASequence throwing an error if the sequence contains characters other
than A, C, G, T, N.
"""
function sequence(bam::BAMAlignment)
    # This follows closely from @encode_seq, but unfortunately not enough to
    # reuse that.
    seq_length = bam.seq_length
    seqdata = zeros(UInt64, Seq.seq_data_len(seq_length))
    ns = falses(Int(seq_length))
    recode_bam_sequence!(bam.data.data, Int(bam.seq_position), Int(bam.seq_length), seqdata, ns)
    return DNASequence(seqdata, ns, 1:seq_length, false, false)
end
=#

# TODO: DNA_Gap
# "=ACMGRSVTWYHKDBN" -> [0,16)
const bam_nucs = [
    DNA_Gap, DNA_A, DNA_C, DNA_M,
    DNA_G,   DNA_R, DNA_S, DNA_V,
    DNA_T,   DNA_W, DNA_Y, DNA_H,
    DNA_K,   DNA_D, DNA_B, DNA_N
]

function sequence(aln::BAMAlignment)
    seq = DNASequence(aln.seq_length)
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

function sequence!(bam::BAMAlignment, seq::DNASequence)
    if !seq.mutable
        error("Cannot copy sequence from BAMAlignment to immutable DNASequence. Call `mutable!(seq)` first.")
    end

    len = Seq.seq_data_len(bam.seq_length)
    if length(seq.data) < len
        resize!(seq.data, len)
    end

    if length(seq.ns) < bam.seq_length
        resize!(seq.ns, bam.seq_length)
    end

    fill!(seq.data, 0)
    fill!(seq.ns, false)
    seq.part = 1:bam.seq_length

    recode_bam_sequence!(bam.data.data, Int(bam.seq_position), Int(bam.seq_length), seq.data, seq.ns)

    return seq
end


# Recode a BAM sequence from 4-bit to 2-bit.
function recode_bam_sequence!(input::Vector{UInt8}, seq_position::Int,
                              seq_length::Int, seqdata::Vector{UInt64}, ns::BitVector)
    len = Seq.seq_data_len(seq_length)
    ored_nucs = UInt8(0)
    j = seq_position
    idx = 1
    for i in 1:len
        shift = 0
        data_i = UInt64(0)
        while shift < 64
            if idx > seq_length
                break
            end

            # high nibble
            @inbounds nt = bam4bit_to_dna[1 + (input[j] >>> 4)]
            if nt == DNA_N
                d = (idx - 1) >>> 6
                r = (idx - 1) & 63
                ns.chunks[d + 1] |= UInt64(1) << r
            else
                ored_nucs |= convert(UInt8, nt)
                data_i |= convert(UInt64, nt) << shift
            end

            shift += 2
            idx += 1

            if idx > seq_length
                break
            end

            # low nibble
            @inbounds nt = bam4bit_to_dna[1 + (input[j] & 0xf)]
            if nt == DNA_N
                d = (idx - 1) >>> 6
                r = (idx - 1) & 63
                ns.chunks[d + 1] |= UInt64(1) << r
            else
                ored_nucs |= convert(UInt8, nt)
                data_i |= convert(UInt64, nt) << shift
            end

            shift += 2
            idx += 1
            j += 1
        end
        @inbounds seqdata[i] = data_i
    end

    # if there was a bad nucleotide, go back and find it
    if ored_nucs & 0b1000 != 0
        idx = 1
        j = seq_position
        while true
            if idx > seq_length
                break
            end

            # high nibble
            nt = bam4bit_to_dna[1 + (input[j] >>> 4)]
            if nt == DNA_INVALID
                error(string(bam4bit_to_char[1 + (input[j] >>> 4)], " is not a valid DNA nucleotide."))
            end
            idx += 1

            if idx > seq_length
                break
            end

            # low nibble
            nt = bam4bit_to_dna[1 + (input[j] & 0xf)]
            if nt == DNA_INVALID
                error(string(bam4bit_to_char[1 + (input[j] & 0xf)], " is not a valid DNA nucleotide."))
            end
            idx += 1
            j += 1
        end
    end
end


function qualities(bam::BAMAlignment)
    r = bam.qual_position:bam.aux_position-1
    qs = Vector{Int8}(length(r))
    for (i, j) in enumerate(r)
        @inbounds qs[i] = bam.data[j] - 33
    end
    return qs
end


function qualities!(bam::BAMAlignment, qs::Vector{Int8})
    r = bam.qual_position:bam.aux_position-1
    if length(qs) != length(r)
        resize!(qs, length(r))
    end
    for (i, j) in enumerate(r)
        @inbounds qs[i] = bam.data.data[j] - 33
    end
    return qs
end


# Return all auxillary data
function auxiliary(bam::BAMAlignment)
    aux = Dict{Tuple{Char, Char}, Any}()
    data = bam.data.data
    p = Int(bam.aux_position)
    while p <= bam.data.part.stop
        tag = (Char(data[p]), Char(data[p + 1]))
        p += 2
        if data[p] == 'B'
            typ = Vector{auxtype[data[p+1]]}
            p += 2
        else
            typ = auxtype[data[p]]
            p += 1
        end
        p, value = getauxdata(data, p, typ)
        aux[tag] = value
    end

    return aux
end


const auxtype = Dict{Char, Type}(
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


# Return a specific tag
function auxillary(bam::BAMAlignment, t1::Char, t2::Char)
    data = bam.data.data
    p = Int(bam.aux_position)

    while p <= bam.data.part.stop && (data[p] != t1 || data[p + 1] != t2)
        p = next_tag_position(data, p)
    end

    if p > bam.data.part.stop
        throw(KeyError(string(t1, t2)))
    else
        p += 2
        if data[p] == 'B'
            typ = Vector{auxtype[data[p+1]]}
            p += 2
        else
            typ = auxtype[data[p]]
            p += 1
        end
        p, value = getauxdata(data, p, typ)
        return value
    end
end


function getauxdata{T}(data::Vector{UInt8}, p::Int, ::Type{T})
    return (Int(p + sizeof(T)), unsafe_load(Ptr{T}(pointer(data, p))))
end


function getauxdata(data::Vector{UInt8}, p::Int, ::Type{Char})
    return (p + 1, Char(unsafe_load(pointer(data, p))))
end


function getauxdata{T}(data::Vector{UInt8}, p::Int, ::Type{Vector{T}})
    n = unsafe_load(Ptr{Int32}(pointer(data, p)))
    p += 4
    xs = Array(T, n)
    unsafe_copy!(pointer(xs), Ptr{T}(pointer(data, p)), n)
    return (Int(p + n * sizeof(T)), xs)
end


function getauxdata(data::Vector{UInt8}, p::Int, ::Type{ASCIIString})
    dataptr = pointer(data, p)
    endptr = ccall(:memchr, Ptr{Void}, (Ptr{Void}, Cint, Csize_t),
                   dataptr, '\0', length(data) - p + 1)
    q = p + (endptr - dataptr) - 1
    return (Int(q + 2), ASCIIString(data[p:q]))
end


function next_tag_position(data::Vector{UInt8}, p::Int)
    typ = data[p + 2]
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
        while data[p] != '\0'
            p += 1
        end
        p += 1
    elseif typ == 'B'
        eltyp = data[p]
        elsize = eltyp == 'c' || eltyp == 'C' ? 1 :
                 eltyp == 's' || eltyp == 'S' ? 2 :
                 eltyp == 'i' || eltye == 'I' || eltyp == 'f' ? 4 :
                 error(string("Unrecognized auxiliary type ", eltyp))
        p += 1
        n = unsafe_load(Ptr{Int32}(pointer(data, p)))
        p += elsize * n
    else
        error(string("Unrecognized auxiliary type ", typ))
    end
    return p
end


function auxiliary(bam::BAMAlignment, tag::AbstractString)
    if length(tag) != 2
        error("BAM auxillary data tags must be of length 2")
    end
    return auxillary(bam, tag[1], tag[2])
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


# FIXME: <: AbstractParser{BAMAlignment}
immutable BAMParser{T <: BufferedInputStream} <: Bio.AbstractParser
    stream::T
    header_text::StringField
    refs::Vector{Tuple{StringField, Int}}
end


@inline function unsafe_load_type(buffer::Vector{UInt8}, position::Integer, T::Type)
    return unsafe_load(convert(Ptr{T}, pointer(buffer, position)))
end


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


function Base.read!(parser::BAMParser, aln::BAMAlignment)
    stream = parser.stream
    if eof(stream)
        # FIXME: throw EOFError
        return false
    end

    # FIXME: why "read the entire entry"?
    # read the entire entry into the buffer
    block_size = read(stream, Int32)
    anchor!(stream)
    seekforward(stream, block_size)
    p = upanchor!(stream)

    # extract fields
    ptr = pointer(stream.buffer, p)
    # FIXME: this may be important in terms of performance: read fixed-length
    # fields into memory
    fields = unsafe_load(convert(Ptr{BAMEntryHead}, ptr))
    @inbounds if 0 <= fields.refid < length(parser.refs)
        aln.seqname = parser.refs[fields.refid + 1][1]
    else
        empty!(aln.seqname)
    end

    aln.position =  fields.pos + 1  # make 1-based
    l_read_name  =  fields.bin_mq_nl        & 0xff
    aln.mapq     = (fields.bin_mq_nl >>  8) & 0xff
    aln.bin      =  fields.bin_mq_nl >> 16
    n_cigar_op   =  fields.flag_nc          & 0xffff
    aln.flag     =  fields.flag_nc   >> 16

    @inbounds if 0 <= fields.next_refid < length(parser.refs)
        aln.next_seqname = parser.refs[fields.next_refid + 1][1]
    else
        empty!(aln.next_seqname)
    end

    aln.next_pos = fields.next_pos
    aln.tlen = fields.tlen

    p += sizeof(BAMEntryHead)
    datalen = stream.position - p
    if length(aln.data) < datalen
        resize!(aln.data, datalen)
    end
    copy!(aln.data, 1, stream.buffer, p, datalen)

    aln.cigar_position = l_read_name + 1
    aln.seq_position   = aln.cigar_position + 4n_cigar_op
    aln.seq_length     = fields.l_seq
    aln.qual_position  = aln.seq_position + cld(fields.l_seq, 2)
    aln.aux_position   = aln.qual_position + fields.l_seq

    return true
end
