# Bio.IO
# ======
#
# I/O interfaces.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module IO

"""
Abstract file format type.

See `subtypes(FileFormat)` for all available file formats.
"""
abstract FileFormat

"""
Abstract formatted input/output type.
"""
abstract AbstractFormattedIO

"""
    stream(io::AbstractFormattedIO)

Return the underlying `IO` object; subtypes of `AbstractFormattedIO` must
implement this method.
"""
function stream end

# delegate method call
for f in (:eof, :flush, :close)
    @eval function Base.$(f)(io::AbstractFormattedIO)
        return $(f)(stream(io))
    end
end

"""
Abstract data reader type.

See `subtypes(AbstractReader)` for all available data readers.
"""
abstract AbstractReader <: AbstractFormattedIO

Base.iteratorsize(::AbstractReader) = Base.SizeUnknown()

"""
Abstract data writer type.

See `subtypes(AbstractWriter)` for all available data writers.
"""
abstract AbstractWriter <: AbstractFormattedIO

end  # module Bio.IO
