__precompile__()

module Bio

include("IO.jl")
include("Ragel.jl")
include("StringFields.jl")
include("util/Util.jl")
include("seq/Seq.jl")
include("intervals/Intervals.jl")
include("align/Align.jl")
include("util/tokenize.jl")
include("util/indexing.jl")
include("phylo/Phylo.jl")
include("structure/Structure.jl")
include("tools/Tools.jl")
include("util/windows.jl")
include("precompile.jl")

end  # module Bio
