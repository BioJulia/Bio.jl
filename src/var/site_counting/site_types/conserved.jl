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

# Methods for the naive framework.
# --------------------------------

ischange(::Type{Conserved}, a::BioSequence, b::BioSequence, idx) = issite(Match, a, b, idx)

# Methods for the bitparallel framework.
# --------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline count_bitpar(::Type{Conserved}, ::Type{$A{2}}, a::UInt64, b::UInt64) = count_bitpar(Match, $A{2}, a, b)
        @inline correct_emptyspace(::Type{Conserved}, ::Type{$A{2}}) = true

        @inline function count_bitpar(::Type{Conserved}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            k, c = count_bitpar(Mutated, $A{4}, a, b)
            return c - k, c
        end
    end
end



"""
    create_nibble_mask(::Type{Conserved}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides are the same.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Conserved}, a::UInt64, b::UInt64)
    certainmask = create_nibble_mask(Certain, a, b)
    return create_nibble_mask(a $ b, 0x0000000000000000) & certainmask
end

"""
    count_nibbles(::Type{Conserved}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
Conserved between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.
**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
Conserved. For example, 'A' and 'R', or 'A' and '-' will not be counted.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Conserved}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Conserved, a, b))
end
