__precompile__()

module Bio

import BioSequences
import GenomicFeatures
import BioAlignments
import Phylogenies
import BioStructures
import GeneticVariation
import BioServices

const Seq = BioSequences
const Intervals = GenomicFeatures
const Align = BioAlignments
const Phylo = Phylogenies
const Structure = BioStructures
const Var = GeneticVariation
const Services = BioServices

include("tools/Tools.jl")

end  # module Bio
