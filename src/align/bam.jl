
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
    # TODO
    #aux::
end


function BAMAlignment()
    return BAMAlignment(StringField(), -1, 0, 0, 0, StringField(), -1, -1,
                        StringField(), StringField(), DNASequence(), UInt8[])
end


immutable BAMParser{T <: BufferedInputStream} <: AbstractParser
    stream::T
    header_text::StringField
    refs::Vector{Tuple{StringField, Int}}

    # TODO: Store header somewhere?
end


@inline function unsafe_load_type(buffer::Vector{UInt8}, position::Integer, T::Type)
    return unsafe_load(convert(Ptr{T}, pointer(buffer, position)))
end


function Base.open(source_stream::BufferedInputStream, ::Type{BAM})
    #bgzf_source = BGZFSource(source_stream)
    #stream = BufferedInputStream(bgzf_source)
    stream = BufferedInputStream(BGZFSource(source_stream))
    #stream = ZlibInflateInputStream(source_stream)

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

    return BAMParser(stream, header_text, refs)
end


#function Base.read!(parser::BAMParser, align::BAMAlignment)
    #stream = parser.stream # 0.05
    #if eof(stream)
        #return false
    #end

    ##block_size = read(stream, Int32) # 0.25

    #anchor!(stream)
    #seekforward(stream, 4)
    #block_size = unsafe_load_type(stream.buffer, upanchor!(stream), Int32)

    ### read the entire entry into the buffer
    #anchor!(stream)
    #seekforward(stream, block_size) # 0.62
    #p = upanchor!(stream)

    #return true
#end


# Maybe if we read a bunch of things at once?
function Base.read!(parser::BAMParser, align::BAMAlignment)
    stream = parser.stream
    while !eof(stream)
        anchor!(stream)
        seekforward(stream, 4)
        block_size = unsafe_load_type(stream.buffer, upanchor!(stream), Int32)

        anchor!(stream)
        seekforward(stream, block_size)
        p = upanchor!(stream)
    end

    return false
end


