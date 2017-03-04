# VCF Writer
# ==========

"""
    VCFWriter(output::IO, header::VCFHeader)

Create a data writer of the VCF file format.

# Arguments
* `output`: data sink
* `header`: VCF header object
"""
type VCFWriter{T<:IO} <: Bio.IO.AbstractWriter
    stream::T
end

function VCFWriter(output::IO, header::VCFHeader)
    writer = VCFWriter(output)
    write(writer, header)
    return writer
end

function Bio.IO.stream(writer::VCFWriter)
    return writer.stream
end

function Base.write(writer::VCFWriter, header::VCFHeader)
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

function Base.write(writer::VCFWriter, record::VCFRecord)
    return write(writer.stream, record, '\n')
end

