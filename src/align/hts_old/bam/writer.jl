# BAM Writer
# ==========

"""
    BAMWriter(output::BGZFStream, header::SAMHeader)

Create a data writer of the BAM file format.

# Arguments
* `output`: data sink
* `header`: SAM header object
"""
type BAMWriter <: Bio.IO.AbstractWriter
    stream::BGZFStreams.BGZFStream
end

function BAMWriter(stream::BGZFStreams.BGZFStream, header::SAMHeader)
    refseqnames = String[]
    refseqlens = Int[]
    for rec in header["SQ"]
        push!(refseqnames, rec["SN"])
        len = rec["LN"]
        if isa(len, AbstractString)
            push!(refseqlens, parse(Int, len))
        else
            push!(refseqlens, len)
        end
    end
    write_header(stream, header, refseqnames, refseqlens)
    return BAMWriter(stream)
end

function Bio.IO.stream(writer::BAMWriter)
    return writer.stream
end

function write_header(stream, header, refseqnames, refseqlens)
    @assert length(refseqnames) == length(refseqlens)
    n = 0

    # magic bytes
    n += write(stream, "BAM\1")

    # SAM header
    buf = IOBuffer()
    l = write(SAMWriter(buf), header)
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

function Base.write(writer::BAMWriter, record::BAMRecord)
    n = 0
    n += unsafe_write(writer.stream, pointer_from_objref(record), BAM_FIXED_FIELDS_BYTES)
    n += unsafe_write(writer.stream, pointer(record.data), data_size(record))
    return n
end
