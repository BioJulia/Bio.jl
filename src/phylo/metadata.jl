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

immutable BasicBranch{T<:AbstractFloat} <: BranchMetaData
    len::Nullable{T}
    conf::Nullable{T}
end

@inline function BasicBranch{T<:AbstractFloat}()
    return BasicBranch(Nullable{T}(), Nullable{T}())
end

@inline function BasicBranch{T<:AbstractFloat}(x::BasicBranch{T}, len::BLWrapper)
    return BasicBranch{B}(convert(Nullable{T}, len.x), x.conf)
end

function branchlength(x::BasicBranch)
    return x.len
end
