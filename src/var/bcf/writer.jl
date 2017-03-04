# BCF Writer
# ==========

"""
    BCFWriter(output::IO, header::VCFHeader)

Create a data writer of the BCF file format.

# Arguments
* `output`: data sink
* `header`: VCF header object
"""
type BCFWriter{T<:IO} <: Bio.IO.AbstractWriter
    stream::BGZFStream{T}
end

function BCFWriter(output::IO, header::VCFHeader)
    stream = BGZFStream(output, "w")
    write(stream, b"BCF\x02\x02")
    buf = IOBuffer()
    len = write(buf, header)
    if len > typemax(Int32)
        error("too long header")
    end
    write(stream, htol(Int32(len)))
    data = takebuf_array(buf)
    @assert length(data) == len
    write(stream, data)
    return BCFWriter(stream)
end

function Bio.IO.stream(writer::BCFWriter)
    return writer.stream
end

function Base.write(writer::BCFWriter, record::BCFRecord)
    n = 0
    n += write(writer.stream, htol(record.sharedlen))
    n += write(writer.stream, htol(record.indivlen))
    n += write(writer.stream, record.data)
    return n
end
