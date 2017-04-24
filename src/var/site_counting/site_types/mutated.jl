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
            m = create_nibble_mask(Certain, a, b)
            d = (a $ b) & m
            return count_nonzero_nibbles(d), count_one_nibbles(m)
        end
    end
end

"""
    count_nibbles(::Type{Mutated}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
any mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.
**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Mutated}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Mutated, a, b))
end

"""
    create_nibble_mask(::Type{Mutated}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides are different.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Mutated}, a::UInt64, b::UInt64)
    certainmask = create_nibble_mask(Certain, a, b)
    return (~create_nibble_mask(a $ b, 0x0000000000000000)) & certainmask
end
