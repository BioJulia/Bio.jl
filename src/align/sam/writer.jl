# SAM Writer
# ==========

"""
    Writer(output::IO, header::Header=Header())

Create a data writer of the SAM file format.

# Arguments
* `output`: data sink
* `header=Header()`: SAM header object
"""
type Writer <: Bio.IO.AbstractWriter
    stream::IO

    function Writer(output::IO, header::Header=Header())
        writer = new(output)
        write(writer, header)
        return writer
    end
end

function Bio.IO.stream(writer::Writer)
    return writer.stream
end

function Base.write(writer::Writer, header::Header)
    n = 0
    for metainfo in header
        n += write(writer, metainfo)
    end
    return n
end

function Base.write(writer::Writer, metainfo::MetaInfo)
    checkfilled(metainfo)
    return write(writer.stream, metainfo, '\n')
end

function Base.write(writer::Writer, record::Record)
    checkfilled(record)
    return write(writer.stream, record, '\n')
end
