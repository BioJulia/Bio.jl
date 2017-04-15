# certain.jl
# ==========
#
# Define Certain sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol or a gap.
"""
immutable Certain <: Site end
const CERTAIN = Certain()

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleic acid in a sequence is certain."
issite{T<:NucleicAcid}(::Type{Certain}, a::T) = iscertain(a)

"Test whether a nucleotide site of two aligned sequences has two certain NucleicAcids."
@inline function issite{T<:NucleicAcid}(::Type{Certain}, a::T, b::T)
    return issite(Certain, a) & issite(Certain, b)
end

# Methods for the bitparallel framework.
# --------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function count_algorithm(s::Certain, a::BioSequence{$A{2}}, b::BioSequence{$A{2}})
            return ALL
        end

        @inline function count_algorithm(s::Certain, a::BioSequence{$A{4}}, b::BioSequence{$A{4}})
            return BITPAR
        end

        @inline function count_bitpar(::Type{Certain}, ::Type{$A{4}}, x::UInt64)
            x = enumerate_nibbles(x)
            x $= 0x1111111111111111
            return count_zero_nibbles(x)
        end

        @inline function count_bitpar(::Type{Certain}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            x = enumerate_nibbles(a) $ 0x1111111111111111
            y = enumerate_nibbles(b) $ 0x1111111111111111
            return count_zero_nibbles(x | y)
        end
    end
end

"""
    create_nibble_mask(::Type{Certain}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites that should be
considered, when counting pairwise mutations between sequences.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(::Type{Certain}, x::UInt64)
    return nibble_mask(enumerate_nibbles(x), 0x1111111111111111)
end

"""
    nibble_mask(::Type{Certain}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites that should be
considered during pairwise distance computation at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(::Type{Certain}, a::UInt64, b::UInt64)
    return nibble_mask(Certain, a) & nibble_mask(Certain, b)
end
