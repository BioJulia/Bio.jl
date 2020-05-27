__precompile__()

module Bio

import BioSequences
import GenomicFeatures
import BioAlignments
import BioStructures
import GeneticVariation
import BioServices
import BioTools

# TODO
#import Phylogenies

const Seq = BioSequences
const Intervals = GenomicFeatures
const Align = BioAlignments
const Structure = BioStructures
const Var = GeneticVariation
const Services = BioServices
const Tools = BioTools

# TODO
#const Phylo = Phylogenies

function __init__()
    @warn "This package has been depreceated, and should not be used for new projects. Please see Bio.jl's repository README for more information."
end

end  # module Bio
