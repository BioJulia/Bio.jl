# mutated.jl
# ============
#
# Define mutated sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Mutated` site describes a site where two aligned nucleotides are definately
mutated. By definately mutated this means that the symbols of the site are
non-ambiguity symbols, and they are not the same symbol.
"""
immutable Mutated <: Mutation end

# Methods for the naive framework.
# --------------------------------

@inline function ischange{T<:NucleicAcid}(::Type{Mutated}, a::T, b::T)
    return issite(Mismatch, a, b)
end

# Methods for the bitparallel framework.
# --------------------------------------
