module Bio

abstract AbstractParser
abstract FileFormat

include("StringFields.jl")
include("Ragel.jl")
include("seq/Seq.jl")
include("intervals/Intervals.jl")
include("align/Align.jl")
include("var/Var.jl")

end # module Bio
