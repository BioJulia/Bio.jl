module Bio

abstract AbstractParser
abstract FileFormat

include("BGZF.jl")
include("StringFields.jl")
include("Ragel.jl")
include("seq/Seq.jl")
include("align/Align.jl")
include("intervals/Intervals.jl")

end # module Bio
