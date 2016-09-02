# Writer
# ======

"""
    FASTQWriter(output::IO; quality_header=false, quality_encoding=:sanger)

Create a data writer of the FASTQ file format.

# Arguments
* `output`: data sink
* `quality_header=false`: output the title line at the third line just after '+'
* `quality_encoding=:sanger`: encoding of base qualities; see `FASTQReader` for available values
"""
type FASTQWriter{T<:IO} <: Bio.IO.AbstractWriter
    output::T

    # write sequence name and description header at the third line
    quality_header::Bool

    # quality encoding used for writing
    quality_encoding::QualityEncoding
end

function Bio.IO.stream(writer::FASTQWriter)
    return writer.output
end

function FASTQWriter(output::IO;
                     quality_header::Bool=false,
                     quality_encoding::Symbol=:sanger)
    return FASTQWriter(output, quality_header, sym2qualenc(quality_encoding))
end

function Base.write(writer::FASTQWriter, seqrec::FASTQSeqRecord)
    output = writer.output
    n = 0

    # header
    n += write(output, '@', seqrec.name)
    if !isempty(seqrec.metadata.description)
        n += write(output, ' ', seqrec.metadata.description)
    end
    n += write(output, '\n')

    # sequence
    for x in seqrec.seq
        n += write(output, Char(x))
    end
    n += write(output, '\n')

    # (optional) identifier and description
    n += write(output, '+')
    if writer.quality_header
        n += write(output, seqrec.name)
        if !isempty(seqrec.metadata.description)
            n += write(output, ' ', seqrec.metadata.description)
        end
    end
    n += write(output, '\n')

    for q in seqrec.metadata.quality
        n += write(output, Char(q + ascii_offset(writer.quality_encoding)))
    end
    n += write(output, '\n')

    return n
end
