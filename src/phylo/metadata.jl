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
# should be of type Nullable.

immutable BasicBranch{T<:AbstractFloat} <: BranchMetaData
    len::Nullable{T}
    conf::Nullable{T}
end

@inline function BasicBranch{T<:AbstractFloat}()
    return BasicBranch(Nullable{T}(), Nullable{T}())
end

@inline function new_branchlength{T<:AbstractFloat,N<:AbstractFloat}(x::BasicBranch{T}, len::Nullable{N})
    return BasicBranch(convert(Nullable{T}, len), x.conf)
end

branchlength(x::BasicBranch) = x.len
