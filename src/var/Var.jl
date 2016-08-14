# The Bio.Var module
# ====================
#
# Types and methods for working with genetic variance,
# including SNPs, Genotypes, and other kinds of population genetic data types.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md



module Var

using DataArrays, Bio.Indexers

import Base:
    convert,
    show


end

include("binary_snps.jl")
