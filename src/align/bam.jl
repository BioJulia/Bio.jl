
using Libz
using BufferedStreams

using Bio.Seq
using Bio.BGZF
using Bio.StringFields

immutable BAM <: FileFormat end


# TODO: This is a temporary type. We should find a more general alignment
# representation
type BAMAlignment
    refname::StringField
    position::Int64
    mapq::UInt8
    bin::UInt16
    flag::UInt16
    next_refname::StringField
    next_pos::Int64
    tlen::Int32
    read_name::StringField
    cigar::StringField # TODO: a more specialized structure
    seq::DNASequence
    qual::Vector{UInt8}
    aux::StringField # TODO: a more specialized structure
end


function BAMAlignment()
    return BAMAlignment(StringField(), -1, 0, 0, 0, StringField(), -1, -1,
                        StringField(), StringField(), DNASequence(), UInt8[],
                        StringField())
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
    if fields.refid >= 0
        alignment.refname = parser.refs[fields.refid + 1][1]
    else
        empty!(alignment.refname)
    end

    alignment.position = fields.pos + 1 # make 1-based
    l_read_name = fields.bin_mq_nl & 0xf
    alignment.mapq = (fields.bin_mq_nl >> 8) & 0xf
    alignment.bin = (fields.bin_mq_nl >> 16) & 0xff
    n_cigar_op = fields.flag_nc & 0xff
    alignment.flag = fields.flag_nc >> 16

    if fields.next_refid >= 0
        alignment.next_refname = parser.refs[fields.next_refid + 1][1]
    else
        empty!(alignment.next_refname)
    end

    alignment.next_pos = fields.next_pos
    alignment.tlen = fields.tlen

    p += sizeof(BAMEntryHead)
    copy!(alignment.read_name, stream.buffer, p, p + l_read_name - 1)

    p += l_read_name
    copy!(alignment.cigar, stream.buffer, p, p + n_cigar_op * 4 - 1)

    # TODO: BAM 4-bit encodes a DNA sequence. We need to write a special
    # conversion function.
    p += n_cigar_op * 4
    #copy!(alignment.seq, stream.buffer, p, p + fields.l_seq - 1)

    p += div(fields.l_seq + 1, 2)
    if length(alignment.qual) != fields.l_seq
        resize!(alignment.qual, fields.l_seq)
    end
    copy!(alignment.qual, 1, stream.buffer, p, fields.l_seq)

    p += fields.l_seq
    copy!(alignment.aux, stream.buffer, p, stream.position - 1)

    return true
end



