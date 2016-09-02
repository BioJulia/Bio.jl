# Writer
# ======

type FASTQWriter{T<:IO} <: Bio.IO.AbstractWriter
    output::T

    # write sequence name and description header at the third line
    quality_header::Bool

    # encoding offset of base qualities (ASCII code = Phred quality + offset)
    ascii_offset::Int
end

function Bio.IO.stream(writer::FASTQWriter)
    return writer.output
end

function FASTQWriter(output::IO, quality_encoding::QualityEncoding;
                     quality_header::Bool=false)
    offset = ascii_encoding_offsets[quality_encoding]
    return FASTQWriter(output, quality_header, offset)
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
        n += write(output, Char(q + writer.ascii_offset))
    end
    n += write(output, '\n')

    return n
end

function Base.show(io::IO, seqrec::FASTQSeqRecord)
    write(io, "@", seqrec.name, " ", seqrec.metadata.description, "\n")
    for c in seqrec.seq
        print(io, c)
    end
    write(io, '\n')
    # print quality scores as a unicode bar chart
    for q in seqrec.metadata.quality
        if q <= 0
            write(io, '▁')
        elseif q <= 6
            write(io, '▂')
        elseif q <= 12
            write(io, '▃')
        elseif q <= 18
            write(io, '▄')
        elseif q <= 24
            write(io, '▅')
        elseif q <= 30
            write(io, '▆')
        elseif q <= 36
            write(io, '▇')
        else
            write(io, '█')
        end
    end
    write(io, '\n')
    return
end

function Base.write(io::IO, seqrec::FASTQSeqRecord;
                    offset::Integer=-1, qualheader::Bool=false)

    # choose offset automatically
    if offset < 0
        if !isempty(seqrec.metadata.quality) && minimum(seqrec.metadata.quality) < 0
            offset = 64 # solexa quality offset
        else
            offset = 33  # sanger
        end
    end

    write(io, "@", seqrec.name)
    if !isempty(seqrec.metadata.description)
        write(io, " ", seqrec.metadata.description)
    end
    write(io, "\n")

    for c in seqrec.seq
        print(io, c)
    end
    write(io, "\n")

    write(io, "+")
    if qualheader
        write(io, seqrec.name)
        if !isempty(seqrec.metadata.description)
            write(io, " ", seqrec.metadata.description)
        end
    end
    write(io, "\n")

    for q in seqrec.metadata.quality
        write(io, Char(q + offset))
    end
    write(io, "\n")
end
