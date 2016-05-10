__precompile__()

module Bio

abstract AbstractParser
abstract FileFormat

include("StringFields.jl")
include("Ragel.jl")
include("seq/Seq.jl")
include("intervals/Intervals.jl")
include("align/Align.jl")
include("structure/Structure.jl")
include("tools/blast/Blast.jl")
include("util/indexing.jl")
include("util/windows.jl")
include("precompile.jl")

end # module Bio
