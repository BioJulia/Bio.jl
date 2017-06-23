__precompile__()

module Bio

include("declare.jl")
include("Exceptions.jl")
include("IO.jl")
include("Mem.jl")
include("Ragel.jl")
include("ReaderHelper.jl")
include("RecordHelper.jl")
include("StringFields.jl")
include("util/Util.jl")
include("Seq.jl")
include("services/Services.jl")
include("intervals/Intervals.jl")
include("align/Align.jl")
include("util/tokenize.jl")
include("util/indexing.jl")
include("util/windows.jl")
include("var/Var.jl")
include("Phylo.jl")
include("structure/Structure.jl")
include("tools/Tools.jl")

end  # module Bio
