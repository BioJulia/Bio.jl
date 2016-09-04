# FASTA Writer
# ============

"""
    FASTAWriter(output::IO; width=70)

Create a data writer of the FASTA file format.

# Arguments
* `output`: data sink
* `width=70`: wrapping width of sequence characters
"""
type FASTAWriter{T<:IO} <: Bio.IO.AbstractWriter
    output::T
    # maximum sequence width (no limit when width ≤ 0)
    width::Int
end

function Bio.IO.stream(writer::FASTAWriter)
    return writer.output
end

function FASTAWriter(output::IO; width::Integer=70)
    return FASTAWriter(output, width)
end

function Base.write(writer::FASTAWriter, seqrec::FASTASeqRecord)
    output = writer.output
    n = 0

    # header
    n += write(output, '>', seqrec.name)
    if !isempty(seqrec.metadata.description)
        n += write(output, ' ', seqrec.metadata.description)
    end
    n += write(output, '\n')

    # sequence
    w = writer.width
    for x in seqrec.seq
        if writer.width > 0 && w == 0
            n += write(output, '\n')
            w = writer.width
        end
        n += write(output, Char(x))
        w -= 1
    end
    n += write(output, '\n')

    return n
end
