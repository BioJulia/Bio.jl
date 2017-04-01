# The Bio.Phylo module
# ====================
#
# A Bio.jl module re exporting the Phylogenies package of the BioJulia
# ecosystem.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Phylo

using Reexport

@reexport using Phylogenies

end # module Phylo
