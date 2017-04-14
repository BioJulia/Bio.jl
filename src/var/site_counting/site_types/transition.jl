# transition.jl
# =============
#
# Define transition sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Transition` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transition mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transition mutation: i.e. A<->G, or C<->T.
"""
immutable Transition <: Mutation end
const TRANSITION = Transition()

# Methods for the naive framework.
# --------------------------------

@inline function ischange{T<:NucleicAcid}(::Type{Transition}, a::T, b::T)
    return issite(Mismatch, a, b) & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

# Methods for the bitparallel framework.
# --------------------------------------
