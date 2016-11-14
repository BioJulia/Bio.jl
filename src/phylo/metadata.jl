# phylo/metadata.jl
# =================
#
# Types and methods for metadata attached to phylogenies.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Any type used as branch metadata should have the following methods defined.

# An empty constructor.
# A constructor that takes a BranchLength type as a second argument.

abstract BranchMetaData

immutable BasicBranch <: BranchMetaData
    len::Nullable{Float64}
    conf::Nullable{Float64}
end

@inline function BasicBranch()
    return BasicBranch(Nullable{Float64}(), Nullable{Float64}())
end

@inline function BasicBranch{T}(x::BasicBranch, len::BranchLength{T})
    return BasicBranch(convert(Nullable{Float64}, len), x.conf)
end

function branchlength(x::BasicBranch)
    return x.len
end
