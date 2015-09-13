module Bio

abstract AbstractParser
abstract FileFormat

include("stringfields.jl")
include("ragel.jl")
include("seq/seq.jl")
include("intervals/intervals.jl")

end # module Bio
