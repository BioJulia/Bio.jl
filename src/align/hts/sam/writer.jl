# SAM Writer
# ==========

"""
    SAMWriter(output::IO, header::SAMHeader=SAMHeader())

Create a data writer of the SAM file format.

# Arguments
* `output`: data sink
* `header=SAMHeader()`: SAM header object
"""
type SAMWriter{T<:IO} <: Bio.IO.AbstractWriter
    stream::T
end

function Bio.IO.stream(writer::SAMWriter)
    return writer.stream
end

function SAMWriter(output::IO, header::SAMHeader=SAMHeader())
    writer = SAMWriter(output)
    write(writer, header)
    return writer
end

function Base.write(writer::SAMWriter, header::SAMHeader)
    n = 0
    for metainfo in header
        n += write(writer, metainfo)
    end
    return n
end

function Base.write(writer::SAMWriter, metainfo::SAMMetaInfo)
    checkfilled(metainfo)
    return write(writer.stream, metainfo, '\n')
end

function Base.write(writer::SAMWriter, record::SAMRecord)
    checkfilled(record)
    return write(writer.stream, record, '\n')
end
