__precompile__()

module Bio

import BioSequences
import GenomicFeatures
import BioAlignments
import Phylogenies

const Seq = BioSequences
const Intervals = GenomicFeatures
const Align = BioAlignments
const Phylo = Phylogenies

include("declare.jl")
include("Exceptions.jl")
include("IO.jl")
include("Mem.jl")
include("Ragel.jl")
include("ReaderHelper.jl")
include("RecordHelper.jl")
include("StringFields.jl")
include("util/Util.jl")
include("services/Services.jl")
include("util/tokenize.jl")
include("util/indexing.jl")
include("util/windows.jl")
include("var/Var.jl")
include("structure/Structure.jl")
include("tools/Tools.jl")

end  # module Bio
