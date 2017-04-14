# FASTQ Writer
# ============

immutable Writer <: Bio.IO.AbstractWriter
    output::IO
    quality_header::Bool
end

function Bio.IO.stream(writer::Writer)
    return writer.output
end

"""
    FASTQ.Writer(output::IO; quality_header=false)

Create a data writer of the FASTQ file format.

# Arguments
* `output`: data sink
* `quality_header=false`: output the title line at the third line just after '+'
"""
function Writer(output::IO; quality_header::Bool=false)
    return Writer(output, quality_header)
end

function Base.write(writer::Writer, record::Record)
    checkfilled(record)
    output = writer.output
    n = 0
    # sequence
    n += write(output, '@')
    n += unsafe_write(output, pointer(record.data, first(record.identifier)), length(record.identifier))
    if hasdescription(record)
        n += write(output, ' ')
        n += unsafe_write(output, pointer(record.data, first(record.description)), length(record.description))
    end
    n += write(output, '\n')
    n += unsafe_write(output, pointer(record.data, first(record.sequence)), length(record.sequence))
    n += write(output, '\n')
    # quality
    n += write(output, '+')
    if writer.quality_header
        n += unsafe_write(output, pointer(record.data, first(record.identifier)), length(record.identifier))
        if hasdescription(record)
            n += write(output, ' ')
            n += unsafe_write(output, pointer(record.data, first(record.description)), length(record.description))
        end
    end
    n += write(output, '\n')
    n += unsafe_write(output, pointer(record.data, first(record.quality)), length(record.quality))
    n += write(output, '\n')
    return n
end
