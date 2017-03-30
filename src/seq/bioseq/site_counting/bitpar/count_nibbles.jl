# Bitwise operations for counting sites
# =====================================
#
# Internal bitwise methods used for identifying and counting sites from the
# 2 and 4 bit encoded representations of DNA and RNA sequences.
#
# These functions probably should not be exported, as it's easy to do things
# very wrong with them if you are not careful.
# However, they allow the counting of sites in dna sequences to be wicked fast,
# whilst taking advantage of their compressed representation and avoiding memory
# expensive conversions to large matrices of NucleicAcid types, which take a
# byte of memory each!
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


# Bitparallel site counting - 4 bit encoding
# ==========================================


## Indel sites

"""
    count_nibbles(::Type{Gap}, x::UInt64)
Count the number of gap sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}}
data.
Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Gap}, x::UInt64)
    return count_zero_nibbles(x)
end

"""
    count_nibbles(::Type{Gap}, a::UInt64, b::UInt64)
Count the number of sites in two aligned chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data which contain gap characters.
Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Gap}, a::UInt64, b::UInt64)
    # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
    return count_nibbles(Gap, a) + count_nibbles(Gap, b) - count_nibbles(Gap, a | b)
end


## Certain sites

"""
    count_nibbles(::Type{Certain}, x::UInt64)
An _internal_ function _not for export_, which will count the number of sites in
a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored in
counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Certain}, x::UInt64)
    return count_one_nibbles(create_nibble_mask(Certain, x))
end

"""
    count_nibbles(::Type{Certain}, x::UInt64)
An _internal_ function _not for export_, which will count the number of sites in
aligned chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored
in counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Certain}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Certain, a, b))
end


## Ambiguous sites

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


## Matching sites

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


## Mismatching sites

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
    return 16 - count_nibbles(Match, a, b)
end
