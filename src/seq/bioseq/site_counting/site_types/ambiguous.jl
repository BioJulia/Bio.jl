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

"""
    create_nibble_mask(::Type{Ambiguous}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent ambiguous sites.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Ambiguous}, x::UInt64)
    return ~create_nibble_mask(Certain, x) & ~create_nibble_mask(Indel, x)
end

"""
    create_nibble_mask(::Type{Ambiguous}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have an ambiguous nucleotide at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Ambiguous}, a::UInt64, b::UInt64)
    return create_nibble_mask(Ambiguous, a) | create_nibble_mask(Ambiguous, b)
end

"""
    count_nibbles(::Type{Ambiguous}, x::UInt64)
Count the number of
ambiguous sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data.
Ambiuous sites are defined as those with more than one bit set.
Note here gap - 0000 - then is not ambiguous, even though it is a candidate for
pairwise deletion.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Ambiguous}, x::UInt64)
    return 16 - count_zero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
end

"""
    count_nibbles(::Type{Ambiguous}, a::UInt64, b::UInt64)
Count the number of sites in two aligned chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data which contain ambiguous characters.
Ambiuous sites are defined as those with more than one bit set.
Note here gap - 0000 - then is not ambiguous, even though it is a candidate for
pairwise deletion.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Ambiguous}, a::UInt64, b::UInt64)
    return 16 - count_zero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
end
