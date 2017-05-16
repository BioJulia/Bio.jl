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

# Methods for the naive framework.
# --------------------------------

for Alph in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function ischange{A<:$Alph,B<:$Alph}(::Type{Transition}, a::BioSequence{A}, b::BioSequence{B}, idx)
            ai = a[idx]
            bi = b[idx]
            return (ai != bi) & ((ispurine(ai) & ispurine(bi)) | (ispyrimidine(ai) & ispyrimidine(bi)))
        end
    end
end

# Methods for the bitparallel framework.
# --------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function count_bitpar(::Type{Transition}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            m = nibble_mask(Certain, a, b)
            d = (a $ b) & m
            remcount = count_zero_nibbles(d) # Count of the identical nucleotides.
            tve = count_one_nibbles(nibble_mask(0x9999999999999999, d)) # Count the 1001 transversion edge case.
            d &= (d >> 1)
            d &= 0x7777777777777777
            tvc = count_ones(d)
            tve += (16 - tvc - remcount)
        end
    end
end
