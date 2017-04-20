# conserved.jl
# ============
#
# Define conserved sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Mismatch` site describes a site where two aligned nucleotides are definately
conserved. By definately conserved this means that the symbols of the site are
non-ambiguity symbols, and they are the same symbol.
"""
immutable Conserved <: Mutation end
const CONSERVED = Conserved()

# Methods for the naive framework.
# --------------------------------

ischange(::Type{Conserved}, a::BioSequence, b::BioSequence, idx) = issite(Match, a, b, idx)

# Methods for the bitparallel framework.
# --------------------------------------
