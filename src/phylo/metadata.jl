# phylo/metadata.jl
# =================
#
# Types and methods for metadata attached to phylogenies.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract BranchMetaData

# Every Metadata type should be immutable, and should have the following
# methods defined:
#
# A no-arg constructor. e.g. BasicBranch()
# A method called new_branchlength, which constructs a new metadata variable,
# from the old metadata variable, and a new branchlength. The new branchlength
# should be of type Nullable{Float64}.
# A method branchlength which returns the branchlength

immutable BasicBranch <: BranchMetaData
    len::Nullable{Float64}
    conf::Nullable{Float64}
end

@inline function BasicBranch()
    return BasicBranch(Nullable{Float64}(), Nullable{Float64}())
end

@inline function new_branchlength(x::BasicBranch, len::Nullable{Float64})
    return BasicBranch(len, x.conf)
end

branchlength(x::BasicBranch) = x.len
