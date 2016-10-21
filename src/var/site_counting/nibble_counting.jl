# Bitwise operations for counting mutations
# =========================================
#
# Internal bitwise methods used for identifying and counting mutations from the
# 2 and 4 bit encoded representations of DNA and RNA sequences.
#
# These functions probably should not be exported, as it's easy to do things
# very wrong with them if you are not careful. However, they allow the
# counting of mutations between sequences to be wicked fast, whilst taking
# advantage of their compressed representation and avoiding memory expensive
# conversions to large matrices of Nucleotide types, which take a byte of space
# each!
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


# Bitparallel site counting - 4 bit encoding
# ==========================================


## Indel sites

"""
    count_sites4(::Type{Indel}, x::UInt64)

Count the number of gap sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}}
data.

Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Indel}, x::UInt64)
    return count_zero_nibbles(x)
end

"""
    count_sites4(::Type{Indel}, a::UInt64, b::UInt64)

Count the number of sites in two aligned chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data which contain gap characters.

Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Indel}, a::UInt64, b::UInt64)
    # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
    return count_sites4(Indel, a) + count_sites4(Indel, b) - count_sites4(Indel, a | b)
end


## Certain sites

"""
    count_sites4(::Type{Certain}, x::UInt64)

An _internal_ function _not for export_, which will count the number of sites in
a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored in
counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Certain}, x::UInt64)
    return count_one_nibbles(create_nibble_mask(Certain, x))
end

"""
    count_sites4(::Type{Certain}, x::UInt64)

An _internal_ function _not for export_, which will count the number of sites in
aligned chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored
in counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Certain}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Certain, a, b))
end


## Ambiguous sites

"""
    count_sites4(::Type{Ambiguous}, x::UInt64)

Count the number of
ambiguous sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data.
Ambiuous sites are defined as those with more than one bit set.
Note here gap - 0000 - then is not ambiguous, even though it is a candidate for
pairwise deletion.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Ambiguous}, x::UInt64)
    return 16 - count_zero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
end

"""
    count_sites4(::Type{Ambiguous}, a::UInt64, b::UInt64)

Count the number of sites in two aligned chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data which contain ambiguous characters.
Ambiuous sites are defined as those with more than one bit set.
Note here gap - 0000 - then is not ambiguous, even though it is a candidate for
pairwise deletion.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Ambiguous}, a::UInt64, b::UInt64)
    return 16 - count_zero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
end


## Matching sites

"""
    count_sites4(::Type{Match}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
matching symbols between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** This function will consider unused chunks of ints not yet used by the
BioSequence as matching gap characters, and so this needs to be taken into
account by calling functions.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Match}, a::UInt64, b::UInt64)
    return count_zero_nibbles(a $ b)
end


## Mismatching sites

"""
    count_sites4(::Type{Mismatch}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
matching symbols between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** This function will consider unused chunks of ints not yet used by the
BioSequence as matching gap characters, and so this needs to be taken into
account by calling functions.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Mismatch}, a::UInt64, b::UInt64)
    return 16 - count_sites4(Match, a, b)
end


## Conserved sites

"""
    count_sites4(::Type{Conserved}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
Conserved between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
Conserved. For example, 'A' and 'R', or 'A' and '-' will not be counted.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Conserved}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Conserved, a, b))
end


## Mutated sites

"""
    count_sites4(::Type{Mutated}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
any mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Mutated}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Mutated, a, b))
end

"""
    count_sites4(::Type{Transition}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
transition mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}}
data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Transition}, a::UInt64, b::UInt64)
    pairdelmask = create_nibble_mask(Certain, a, b)
    a &= pairdelmask
    b &= pairdelmask
    diffs = a $ b
    ts_filtered = diffs $ 0xAAAAAAAAAAAAAAAA
    return count_one_nibbles(ts_filtered) + count_one_nibbles(~ts_filtered)
end

"""
    count_sites4(::Type{Transversion}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
transversion mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}}
data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Transversion}, a::UInt64, b::UInt64)
    pairdelmask = create_nibble_mask(Certain, a, b)
    a &= pairdelmask
    b &= pairdelmask
    diffs = a $ b
    tv_filtered_a = diffs $ 0x3333333333333333
    tv_filtered_b = diffs $ 0x6666666666666666
    return count_one_nibbles(tv_filtered_a) + count_one_nibbles(~tv_filtered_a) +
           count_one_nibbles(tv_filtered_b) + count_one_nibbles(~tv_filtered_b)
end
