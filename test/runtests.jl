
const pkglist = [
    "BioSequences",
    "GenomicFeatures",
    "BioAlignments",
#    "Phylogenies",
    "BioStructures",
    #"GeneticVariation",
    "BioServices"
#    "BioTools"
]

import Pkg
for dep in pkglist
    Pkg.test(dep, coverage = false)
end

