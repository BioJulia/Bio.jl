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

function Base.open{T<:AbstractFormattedIO}(f::Function, ::Type{T}, args...; kwargs...)
    io = open(T, args...; kwargs...)
    try
        f(io)
    finally
        close(io)
    end
end


"""
Abstract data reader type.

See `subtypes(AbstractReader)` for all available data readers.
"""
abstract AbstractReader <: AbstractFormattedIO

Base.iteratorsize(::AbstractReader) = Base.SizeUnknown()

function Base.open{T<:AbstractReader}(::Type{T}, filepath::AbstractString, args...; kwargs...)
    return T(open(filepath), args...; kwargs...)
end


"""
Abstract data writer type.

See `subtypes(AbstractWriter)` for all available data writers.
"""
abstract AbstractWriter <: AbstractFormattedIO

function Base.open{T<:AbstractWriter}(::Type{T}, filepath::AbstractString, args...; kwargs...)
    i = findfirst(kwarg -> kwarg[1] == :append, kwargs)
    if i > 0
        append = kwargs[i][1]
        if !isa(append, Bool)
            throw(ArgumentError("append must be boolean"))
        end
    else
        append = false
    end
    return T(open(filepath, append ? "a" : "w"), args...; kwargs...)
end

end  # module Bio.IO
