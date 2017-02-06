# SAM Reader
# ==========

"""
    SAMReader(input::IO)

Create a data reader of the SAM file format.

# Arguments
* `input`: data source
"""
type SAMReader <: Bio.IO.AbstractReader
    state::Ragel.State
    header::SAMHeader

    function SAMReader(input::BufferedStreams.BufferedInputStream)
        reader = new(Ragel.State(samparser_start, input), SAMHeader())
        while !eof(input) && BufferedStreams.peek(input) == UInt8('@')
            # NOTE: This reads a header line, not the first record.
            reader.state.cs = samparser_start
            @assert read(reader) == SAMRecord()
        end
        return reader
    end
end

function SAMReader(input::IO)
    return SAMReader(BufferedStreams.BufferedInputStream(input))
end

function Bio.IO.stream(reader::SAMReader)
    return reader.state.stream
end

function Base.show(io::IO, reader::SAMReader)
    println(io, summary(reader), ":")
      print(io, "  header keys: ", join(keys(reader.header), ", "))
end

function header(reader::SAMReader)
    return reader.header
end

function Base.eltype(::Type{SAMReader})
    return SAMRecord
end
