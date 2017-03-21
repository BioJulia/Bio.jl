import Base.@deprecate, Base.@deprecate_binding
import Base.depwarn

# Bio.jl v0.4
# -----------

# TODO: any way to deprecate intersect?
#@deprecate Base.intersect(bb::BigBedReader, query::Interval)                       eachoverlap(bb, query)
#@deprecate Base.intersect{T}(a::IntervalCollection{T}, b::Interval)                eachoverlap(a, b)
#@deprecate Base.intersect{S,T}(a::IntervalCollection{S}, b::IntervalCollection{T}) eachoverlap(a, b)
#@deprecate Base.intersect(a::IntervalCollection, b::IntervalStreamOrArray)         eachoverlap(a, b)


# Bio.jl v0.3
# -----------

immutable BED <: Bio.IO.FileFormat end
export BED

immutable BigWig <: Bio.IO.FileFormat end
immutable BigBed <: Bio.IO.FileFormat end
export BigWig, BigBed
