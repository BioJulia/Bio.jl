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

ischange(::Type{Mutated}, a::BioSequence, b::BioSequence, idx) = issite(Mismatch, a, b, idx)


# Methods for the bitparallel framework.
# --------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline count_bitpar(::Type{Mutated}, ::Type{$A{2}}, a::UInt64, b::UInt64) = count_bitpar(Mismatch, $A{2}, a, b)

        @inline function count_bitpar(::Type{Mutated}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            m = nibble_mask(Certain, a, b)
            d = (a $ b) & m
            return count_nonzero_nibbles(d), count_one_nibbles(m)
        end
    end
end
