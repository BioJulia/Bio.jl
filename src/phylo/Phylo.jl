# The Bio.Phylo module
# ====================
#
# Types and methods for phylogenetic trees.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md



module Phylo

using Bio.Tokenize

using LightGraphs

using Bio.Indexers

import Bio.FileFormat

export

    # Types
    Phylogeny,

    # Phylogeny Methods

    ## Roots & Clades
    isrooted,
    isrerootable,
    root,
    isroot,
    clades,
    isclade,

    ## Leaves / Taxa
    leaves,
    isleaf,

    ## Structure manipulation
    hasparent,
    parent,
    haschildren,
    nchildren,
    haschild,
    children,
    ispreterminal,
    issemipreterminal,

    ## Misc
    n_possible_rooted,
    n_possible_unrooted



include("phylogeny.jl")
include("node_basics.jl")
include("branch_basics.jl")
include("manipulation.jl")


end # module Phylo
