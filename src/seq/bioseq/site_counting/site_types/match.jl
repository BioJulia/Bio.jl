# match.jl
# ========
#
# Define Matches for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
immutable Match <: Site end
const MATCH = Match()

@inline function count_algorithm{A}(s::Match, a::BioSequence{A}, b::BioSequence{A})
    return BITPAR
end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleotide site of two aligned sequences, constitutes a match."
issite{T<:NucleicAcid}(::Type{Match}, a::T, b::T) = a == b

# Methods for the bitparallel framework.
# --------------------------------------

@inline correct_endspace(::Type{Match}) = true

"""
    nibble_mask(::Type{Match}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the biological
symbols are the same.

Care should be taken as unused sites, as well as indels may be allowed through
this mask as both are seen as a pair of matching 0000 nibbles.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(::Type{Match}, a::UInt64, b::UInt64)
    return nibble_mask(0x0000000000000000, a $ b)
end

"""
    count_nibbles(::Type{Match}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
matching symbols between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.
**Note:** This function will consider unused chunks of ints not yet used by the
BioSequence as matching gap characters, and so this needs to be taken into
account by calling functions.

**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Match}, a::UInt64, b::UInt64)
    return count_zero_nibbles(a $ b)
end

@inline function count_bitpairs(::Type{Match}, a::UInt64, b::UInt64)
    return count_zero_bitpairs(a $ b)
end
