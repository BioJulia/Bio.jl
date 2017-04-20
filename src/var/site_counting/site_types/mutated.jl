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
const MUTATED = Mutated()

# Methods for the naive framework.
# --------------------------------

ischange(::Type{Mutated}, a::BioSequence, b::BioSequence, idx) = issite(Mismatch, a, b, idx)


# Methods for the bitparallel framework.
# --------------------------------------
