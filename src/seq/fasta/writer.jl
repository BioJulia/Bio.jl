# FASTA Writer
# ============

"""
    FASTA.Writer(output::IO; width=70)

Create a data writer of the FASTA file format.

# Arguments
* `output`: data sink
* `width=70`: wrapping width of sequence characters
"""
type Writer <: Bio.IO.AbstractWriter
    output::IO
    # maximum sequence width (no limit when width ≤ 0)
    width::Int
end

function Bio.IO.stream(writer::Writer)
    return writer.output
end

function Writer(output::IO; width::Integer=70)
    return Writer(output, width)
end

function Base.write(writer::Writer, record::Record)
    checkfilled(record)
    output = writer.output
    n::Int = 0
    if writer.width ≤ 0
        n += write(output, record, '\n')
    else
        headerlen = hasdescription(record) ? last(record.description) : last(record.identifier)
        n += unsafe_write(output, pointer(record.data), headerlen)
        n += write(output, '\n')
        p = pointer(record.data, first(record.sequence))
        p_end = pointer(record.data, last(record.sequence))
        while p ≤ p_end
            w = min(writer.width, p_end - p + 1)
            n += unsafe_write(output, p, w)
            n += write(output, '\n')
            p += w
        end
    end
    return n
end
