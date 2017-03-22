# BCF Writer
# ==========

type Writer{T<:IO} <: Bio.IO.AbstractWriter
    stream::BGZFStreams.BGZFStream{T}
end

"""
    BCF.Writer(output::IO, header::VCF.Header)

Create a data writer of the BCF file format.

# Arguments
* `output`: data sink
* `header`: VCF header object
"""
function Writer(output::IO, header::Bio.Var.VCF.Header)
    stream = BGZFStreams.BGZFStream(output, "w")
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
    return Writer(stream)
end

function Bio.IO.stream(writer::Writer)
    return writer.stream
end

function Base.write(writer::Writer, record::Record)
    n = 0
    n += write(writer.stream, htol(record.sharedlen))
    n += write(writer.stream, htol(record.indivlen))
    n += write(writer.stream, record.data)
    return n
end
