# ambiguous.jl
# ============
#
# Define Ambiguous sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
immutable Ambiguous <: Site end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleic acid in a sequence is ambiguous."
issite{T<:NucleicAcid}(::Type{Ambiguous}, a::T) = isambiguous(a)

"Test whether a nucleotide site of two aligned sequences has ambiguities."
@inline function issite{T<:NucleicAcid}(::Type{Ambiguous}, a::T, b::T)
    return issite(Ambiguous, a) | issite(Ambiguous, b)
end

# Methods for the bitparallel framework.
# --------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function count_algorithm(s::Ambiguous, a::BioSequence{$A{2}}, b::BioSequence{$A{2}})
            return NULL
        end

        @inline function count_algorithm(s::Ambiguous, a::BioSequence{$A{4}}, b::BioSequence{$A{4}})
            return BITPAR
        end

        @inline function count_bitpar(::Type{Ambiguous}, ::Type{$A{4}}, x::UInt64)
            return count_nonzero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
        end

        @inline function count_bitpar(::Type{Ambiguous}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_nonzero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
        end
    end
end

"""
    nibble_mask(::Type{Ambiguous}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent ambiguous sites.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(::Type{Ambiguous}, x::UInt64)
    return ~nibble_mask(Certain, x) & ~nibble_mask(Indel, x)
end

"""
    nibble_mask(::Type{Ambiguous}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have an ambiguous nucleotide at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(::Type{Ambiguous}, a::UInt64, b::UInt64)
    return nibble_mask(Ambiguous, a) | nibble_mask(Ambiguous, b)
end
