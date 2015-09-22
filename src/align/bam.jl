
using Libz
using BufferedStreams

using Bio.Seq
using Bio.Seq: DNA_INVALID
using Bio.BGZF
using Bio.StringFields
using Bio.Intervals
using IntervalTrees

immutable BAM <: FileFormat end

# OK, let's really think about how we want to represent BAM.
#
# Maybe I should just store the second half of the BAM file in one StringField
# buffer, then manifest the other fields as needed. Should they be cached?

# TODO: This is a temporary type. We should find a more general alignment
# representation
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
    data::StringField
    cigar_position::Int32
    seq_position::Int32
    seq_length::Int32
    qual_position::Int32
    aux_position::Int32
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
function readname(bam::BAMAlignment)
    bam.data[1:bam.cigar_position-1]
end


"""
Copy the read's name into `name`.
"""
function readname!(bam::BAMAlignment, name::StringField)
    copy!(name, bam.data, 1, bam.cigar_position-1)
end


function cigar(bam::BAMAlignment)
    # TODO
end


const bam4bit_to_char = [
    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N' ]

const bam4bit_to_dna = [
    DNA_INVALID, DNA_A,       DNA_C,       DNA_INVALID,
    DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
    DNA_T,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
    DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_N ]


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
    qs = Array(Int8, length(r))
    for (i, j) in enumerate(r)
        @inbounds qs[i] = bam.data.data[j] - 33
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
    # TODO
end


# Return a specific tag
function auxillary(bam::BAMAlignment, t1::Char, t2::Char)
    # TODO
end


function auxillary(bam::BAMAlignment, tag::AbstractString)
    if length(tag) != 2
        error("BAM auxillary data tags must be of length 2")
    end
    return auxillary(bam, tag[1], tag[2])
end


function BAMAlignment()
    return BAMAlignment(StringField(), -1, 0, 0, 0, StringField(), -1, -1,
                        StringField(), -1, -1, -1, -1, -1)
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


immutable BAMParser{T <: BufferedInputStream} <: AbstractParser
    stream::T
    header_text::StringField
    refs::Vector{Tuple{StringField, Int}}
end


@inline function unsafe_load_type(buffer::Vector{UInt8}, position::Integer, T::Type)
    return unsafe_load(convert(Ptr{T}, pointer(buffer, position)))
end


function Base.open{T <: BufferedInputStream}(source_stream::T, ::Type{BAM})
    stream = BufferedInputStream(BGZFSource(source_stream))

    buffer = stream.buffer

    # magic bytes
    anchor!(stream)
    seekforward(stream, 4)
    p = upanchor!(stream)
    if !(buffer[p]     == UInt8('B') &&
         buffer[p + 1] == UInt8('A') &&
         buffer[p + 2] == UInt8('M') &&
         buffer[p + 3] == 0x01)
        error("Input was not a valid BAM file. (Nonmatching magic bytes.)")
    end

    # header size
    anchor!(stream)
    seekforward(stream, 4)
    header_size = unsafe_load_type(buffer, upanchor!(stream), Int32)

    # header
    anchor!(stream)
    seekforward(stream, header_size)
    header_text = StringField(takeanchored!(stream))
    upanchor!(stream)

    # number of reference sequences
    anchor!(stream)
    seekforward(stream, 4)
    numrefs = unsafe_load_type(buffer, upanchor!(stream), Int32)

    # reference sequence names and sizes
    refs = Array(Tuple{StringField, Int}, numrefs)
    for i in 1:numrefs
        anchor!(stream)
        seekforward(stream, 4)
        name_length = unsafe_load_type(buffer, upanchor!(stream), Int32)

        anchor!(stream)
        seekforward(stream, name_length)
        name = StringField(takeanchored!(stream))

        anchor!(stream)
        seekforward(stream, 4)
        seq_size = unsafe_load_type(buffer, upanchor!(stream), Int32)
        refs[i] = (name, seq_size)
    end

    return BAMParser{BufferedInputStream{BGZFSource{T}}}(stream, header_text, refs)
end


function Base.read!(parser::BAMParser, alignment::BAMAlignment)
    stream = parser.stream
    if eof(stream)
        return false
    end

    # read the entire entry into the buffer
    block_size = read(stream, Int32)
    anchor!(stream)
    seekforward(stream, block_size)
    p = upanchor!(stream)

    # extract fields
    ptr = pointer(stream.buffer, p)
    fields = unsafe_load(convert(Ptr{BAMEntryHead}, ptr))
    @inbounds if 0 <= fields.refid < length(parser.refs)
        alignment.seqname = parser.refs[fields.refid + 1][1]
    else
        empty!(alignment.seqname)
    end

    alignment.position = fields.pos + 1 # make 1-based
    l_read_name = fields.bin_mq_nl & 0xff
    alignment.mapq = (fields.bin_mq_nl >> 8) & 0xff
    alignment.bin = (fields.bin_mq_nl >> 16) & 0xffff
    n_cigar_op = fields.flag_nc & 0xff
    alignment.flag = fields.flag_nc >> 16

    @inbounds if 0 <= fields.next_refid < length(parser.refs)
        alignment.next_seqname = parser.refs[fields.next_refid + 1][1]
    else
        empty!(alignment.next_seqname)
    end

    alignment.next_pos = fields.next_pos
    alignment.tlen = fields.tlen

    p += sizeof(BAMEntryHead)
    copy!(alignment.data, stream.buffer, p, stream.position - 1)

    alignment.cigar_position = l_read_name + 1
    alignment.seq_position = alignment.cigar_position + n_cigar_op * sizeof(Int32)
    alignment.seq_length = fields.l_seq
    alignment.qual_position = alignment.seq_position + div(fields.l_seq + 1, 2)
    alignment.aux_position = alignment.qual_position + fields.l_seq

    return true
end



