# transversion.jl
# ===============
#
# Define transition sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Transversion` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transversion mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transversion mutation: i.e. A<->C,
A<->T, G<->T, G<->C.
"""
immutable Transversion <: Mutation end

# Methods for the naive framework.
# --------------------------------

@inline function ischange{T<:NucleicAcid}(::Type{Transversion}, a::T, b::T)
    return issite(Mismatch, a, b) & ((ispurine(a) & ispyrimidine(b)) | (ispyrimidine(a) & ispurine(b)))
end

# Methods for the bitparallel framework.
# --------------------------------------
