# VCF Writer
# ==========

type Writer{T<:IO} <: Bio.IO.AbstractWriter
    stream::T
end

"""
    VCF.Writer(output::IO, header::VCF.Header)

Create a data writer of the VCF file format.

# Arguments
* `output`: data sink
* `header`: VCF header object
"""
function Writer(output::IO, header::Header)
    writer = Writer(output)
    write(writer, header)
    return writer
end

function Bio.IO.stream(writer::Writer)
    return writer.stream
end

function Base.write(writer::Writer, header::Header)
    n = 0
    for metainfo in header.metainfo
        n += write(writer.stream, metainfo, '\n')
    end
    n += write(writer.stream, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    if !isempty(header.sampleID)
        n += write(writer.stream, "\tFORMAT")
    end
    for id in header.sampleID
        n += write(writer.stream, '\t', id)
    end
    n += write(writer.stream, '\n')
    return n
end

function Base.write(writer::Writer, record::Record)
    return write(writer.stream, record, '\n')
end
