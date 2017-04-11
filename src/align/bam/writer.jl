# BAM Writer
# ==========

"""
    BAM.Writer(output::BGZFStream, header::SAM.Header)

Create a data writer of the BAM file format.

# Arguments
* `output`: data sink
* `header`: SAM header object
"""
type Writer <: Bio.IO.AbstractWriter
    stream::BGZFStreams.BGZFStream
end

function Writer(stream::BGZFStreams.BGZFStream, header::SAM.Header)
    refseqnames = String[]
    refseqlens = Int[]
    for metainfo in find(header, "SQ")
        push!(refseqnames, metainfo["SN"])
        push!(refseqlens, parse(Int, metainfo["LN"]))
    end
    write_header(stream, header, refseqnames, refseqlens)
    return Writer(stream)
end

function Bio.IO.stream(writer::Writer)
    return writer.stream
end

function Base.write(writer::Writer, record::Record)
    n = 0
    n += unsafe_write(writer.stream, pointer_from_objref(record), FIXED_FIELDS_BYTES)
    n += unsafe_write(writer.stream, pointer(record.data), data_size(record))
    return n
end

function write_header(stream, header, refseqnames, refseqlens)
    @assert length(refseqnames) == length(refseqlens)
    n = 0

    # magic bytes
    n += write(stream, "BAM\1")

    # SAM header
    buf = IOBuffer()
    l = write(SAM.Writer(buf), header)
    n += write(stream, Int32(l))
    n += write(stream, takebuf_array(buf))

    # reference sequences
    n += write(stream, Int32(length(refseqnames)))
    for (seqname, seqlen) in zip(refseqnames, refseqlens)
        namelen = length(seqname)
        n += write(stream, Int32(namelen + 1))
        n += write(stream, seqname, '\0')
        n += write(stream, Int32(seqlen))
    end

    return n
end
