# SAM Reader
# ==========

type SAMReader <: Bio.IO.AbstractReader
    state::Ragel.State
    header::SAMHeader

    function SAMReader(input::BufferedInputStream)
        reader = new(Ragel.State(samparser_start, input), SAMHeader())
        while !eof(input) && BufferedStreams.peek(input) == UInt8('@')
            # NOTE: This reads a header line, not the first record.
            @assert read(reader) == SAMRecord()
        end
        return reader
    end
end

function Bio.IO.stream(reader::SAMReader)
    return reader.state.stream
end

function header(reader::SAMReader)
    return reader.header
end

function Base.eltype(::Type{SAMReader})
    return SAMRecord
end

function Base.open(input::BufferedInputStream, ::Type{SAM})
    return SAMReader(input)
end
