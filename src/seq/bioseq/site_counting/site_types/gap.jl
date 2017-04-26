# gap.jl
# =============
#
# Define Gaps for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
An `Gap` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
immutable Gap <: Site end

# Methods for the naive framework.
# --------------------------------

"Test whether a position in a sequence is a gap character."
issite(::Type{Gap}, a::BioSequence, idx) = isgap(a[idx])

"Test whether a position in two aligned sequences has gap characters."
@inline function issite(::Type{Gap}, a::BioSequence, b::BioSequence, idx)
    return issite(Gap, a, idx) | issite(Gap, b, idx)
end






# Methods for the bitparallel framework.
# --------------------------------------

@inline correct_emptyspace{A<:Alphabet}(::Type{Gap}, ::Type{A}) = true

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function count_algorithm(s::Gap, a::BioSequence{$A{2}}, b::BioSequence{$A{2}})
            return NULL
        end

        @inline function count_algorithm(s::Gap, a::BioSequence{$A{4}}, b::BioSequence{$A{4}})
            return BITPAR
        end

        @inline function count_bitpar(::Type{Gap}, ::Type{$A{4}}, x::UInt64)
            return count_zero_nibbles(x)
        end

        @inline function count_bitpar(::Type{Gap}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
            return count_zero_nibbles(a) + count_zero_nibbles(b) - count_zero_nibbles(a | b)
        end
    end
end

"""
    create_nibble_mask(::Type{Gap}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent gaps.

Care should be taken as unused nibbles, as well as indels will be allowed through
this mask.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(::Type{Gap}, x::UInt64)
    return nibble_mask(x, 0x0000000000000000)
end

"""
    create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have a gap at a given nibble.

Care should be taken as unused nibbles, as well as gaps will be allowed through
this mask.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)
    return nibble_mask(Gap, a) | nibble_mask(Gap, b)
end
