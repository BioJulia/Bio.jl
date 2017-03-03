__precompile__()

module Bio

include("declare.jl")
include("IO.jl")
include("Ragel.jl")
include("ReaderHelper.jl")
include("StringFields.jl")
include("util/Util.jl")
include("seq/Seq.jl")
include("services/Services.jl")
include("intervals/Intervals.jl")
include("align/Align.jl")
include("util/tokenize.jl")
include("util/indexing.jl")
include("var/Var.jl")
include("phylo/Phylo.jl")
include("structure/Structure.jl")
include("tools/Tools.jl")
include("util/windows.jl")
include("precompile.jl")

end  # module Bio
