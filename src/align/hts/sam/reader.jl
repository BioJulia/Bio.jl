# SAM Reader
# ==========

type SAMReader <: Bio.IO.AbstractParser
    state::Ragel.State
    header::SAMHeader

    function SAMReader(input::BufferedInputStream)
        return new(Ragel.State(samparser_start, input), SAMHeader())
    end
end

function header(reader::SAMReader)
    return reader.header
end

function Base.eltype(::Type{SAMReader})
    return SAMRecord
end

function Base.eof(reader::SAMReader)
    return eof(reader.state.stream)
end

function Base.close(reader::SAMReader)
    close(reader.state.stream)
end

function Base.open(input::BufferedInputStream, ::Type{SAM})
    return SAMReader(input)
end
