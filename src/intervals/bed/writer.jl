# BED Writer
# ==========

"""
    BED.Writer(output::IO)

Create a data writer of the BED file format.

# Arguments:
* `output`: data sink
"""
immutable Writer <: Bio.IO.AbstractWriter
    output::IO
end

function Bio.IO.stream(writer::Writer)
    return writer.output
end

function Base.write(writer::Writer, record::Record)
    return write(writer.output, record, '\n')
end
