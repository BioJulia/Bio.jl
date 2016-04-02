module Bio

abstract AbstractParser
abstract FileFormat

include("StringFields.jl")
include("Ragel.jl")
include("seq/Seq.jl")
include("util/Util.jl")
include("intervals/Intervals.jl")
include("align/Align.jl")
include("structure/Structure.jl")
include("annotations.jl")
include("tools/tools.jl")
include("phylo/Phylo.jl")

end
