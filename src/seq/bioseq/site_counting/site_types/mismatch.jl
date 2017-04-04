# mismatch.jl
# ===========
#
# Define mismatches for site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
immutable Mismatch <: Site end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleotide site of two aligned sequences, constitutes a mismatch."
issite{T<:NucleicAcid}(::Type{Mismatch}, a::T, b::T) = a != b

# Methods for the bitparallel framework.
# --------------------------------------

"""
    create_nibble_mask(::Type{Mismatch}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the biological
symbols are not the same.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Mismatch}, a::UInt64, b::UInt64)
    return ~create_nibble_mask(Match, a, b)
end

"""
    count_nibbles(::Type{Mismatch}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
matching symbols between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.
**Note:** This function will consider unused chunks of ints not yet used by the
BioSequence as matching gap characters, and so this needs to be taken into
account by calling functions.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Mismatch}, a::UInt64, b::UInt64)
    return count_nonzero_nibbles(a $ b)
end

@inline function count_bitpairs(::Type{Mismatch}, x::UInt64, y::UInt64)
    return count_nonzero_bitpairs(x $ y)
end
