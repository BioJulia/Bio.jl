immutable BED <: Bio.IO.FileFormat end
immutable BigBed <: Bio.IO.FileFormat end
immutable BigWig <: Bio.IO.FileFormat end

export BED

function Base.open{F<:Union{BED,BigBed}}(::AbstractString, ::Type{F})
    error("open(filepath, format) syntax has been removed. Please use open(reader|writer, filepath) instead.")
end

function Base.open{F<:Union{BED}}(::AbstractString, ::AbstractString, ::Type{F})
    error("open(filepath, mode, format) syntax has been removed. Please use open(reader|writer, filepath) instead.")
end
