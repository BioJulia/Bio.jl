# Bio.IO
# ======
#
# I/O interfaces.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module IO

export
    FileFormat,
    AbstractParser,
    AbstractWriter

"""
Abstract file format type.

See `subtypes(FileFormat)` for all available file formats.
"""
abstract FileFormat

"""
Abstract data reader type.

See `subtypes(AbstractParser)` for all available data parsers.
"""
abstract AbstractParser

stream() = nothing

Base.eof{T<:AbstractParser}(parser::T) = eof(stream(parser))
Base.close{T<:AbstractParser}(parser::T) = close(stream(parser))

Base.iteratorsize(::AbstractParser) = Base.SizeUnknown()

"""
Abstract data writer type.

See `subtypes(AbstractWriter)` for all available data writers.
"""
abstract AbstractWriter

stream(parser::AbstractWriter) = nothing
Base.close{T<:AbstractWriter}(writer::T) = close(stream(writer))
Base.flush{T<:AbstractWriter}(writer::T) = flush(stream(writer))

end  # module Bio.IO
