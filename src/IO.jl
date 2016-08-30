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
    AbstractReader,
    AbstractWriter

"""
Abstract file format type.

See `subtypes(FileFormat)` for all available file formats.
"""
abstract FileFormat

"""
Abstract data reader type.

See `subtypes(AbstractReader)` for all available data readers.
"""
abstract AbstractReader

Base.iteratorsize(::AbstractReader) = Base.SizeUnknown()

"""
Abstract data writer type.

See `subtypes(AbstractWriter)` for all available data writers.
"""
abstract AbstractWriter

end  # module Bio.IO
