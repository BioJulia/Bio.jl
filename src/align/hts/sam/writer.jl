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
    if haskey(header, "HD")
        n += write_headerline(writer.stream, "HD", header["HD"])
    end

    for tag in keys(header)
        if tag == "HD"
            continue
        end
        for h in header[tag]
            n += write_headerline(writer.stream, tag, h)
        end
    end

    return n
end

function write_headerline(io, tag, line)
    n = 0
    n += write(io, '@', tag)
    for (t, val) in line
        checkkeytag(t)
        n += write(io, '\t', t, ':', string(val))
    end
    n += write(io, '\n')
    return n
end

function Base.write(writer::SAMWriter, record::SAMRecord)
    stream = writer.stream

    n = 0
    n += write(
        stream,
        record.name, '\t',
        string(record.flag), '\t',
        record.refname, '\t',
        string(record.pos), '\t',
        string(record.mapq), '\t',
        record.cigar, '\t',
        record.next_refname, '\t',
        string(record.next_pos), '\t',
        string(record.tlen), '\t',
        record.seq, '\t',
        record.qual)
    n += write_optfields(stream, record.optional_fields)
    n += write(stream, '\n')

    return n
end

function write_optfields(io, optfields)
    n = 0
    for (tag, val) in optfields
        checkkeytag(tag)
        n += write(io, '\t', tag, ':', auxtypechar[typeof(val)], ':', string(val))
    end
    return n
end
