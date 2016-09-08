# Bitwise operations for identifying and counting mutations
# =========================================================
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


# Nibble operations
# =================

"""
    enumerate_nibbles(abxor::UInt64)

Count the number of set bits in each nibble of a unsigned 64 bit integer
(With BioSequence 4 bit encoding, each nibble is one RNA or DNA nucleotide).

**This is an internal method and should not be exported.**

E.g. An input of:

0100 0010 0001 0110 1100 1110 1101 1111

Would result in:

0001 0001 0001 0010 0010 0011 0011 0100

This is used to identify different occurances of certain bit patterns.
"""
@inline function enumerate_nibbles(x::UInt64)
    x = x - ((x >> 1) & 0x5555555555555555)
    return (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
end

"""
    count_zero_nibbles(x::UInt64)

Counts the number of nibbles in a UInt64 `x` that have all their bits unset i.e.
all nibbles of 0000.

**This is an internal method and should not be exported.**

E.g. An input of:

0x0F11F111F11111F1

Would give the answer: 1.
"""
@inline function count_zero_nibbles(x::UInt64)
    return 16 - count_ones((x & 0x1111111111111111) |
    (x & 0x2222222222222222) >> 1 |
    (x & 0x4444444444444444) >> 2 |
    (x & 0x8888888888888888) >> 3)
end

"""
    count_one_nibbles(x::UInt64)

Counts the number of nibbles in a UInt64 `x` that have all their bits set i.e.
all nibbles of 1111.

**This is an internal method and should not be exported.**

E.g. An input of:

0x00A00C200010000F

Would give the answer: 1.
"""
@inline function count_one_nibbles(x::UInt64)
    return count_ones((x & 0x1111111111111111) &
    ((x & 0x2222222222222222) >> 1) &
    ((x & 0x4444444444444444) >> 2) &
    ((x & 0x8888888888888888) >> 3))
end

# Nibble masking functions
# ------------------------

"""
    create_nibble_mask(x::UInt64, value::UInt64)

Create a mask for the nibbles (groups of four bits) in a 64 bit integer `x`
that match a given value dictated by the pattern in `value`.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(value::UInt64, x::UInt64)
    # XOR with the desired values. So Conserveding nibbles will be 0000.
    x $= value
    # Horizontally OR the nibbles.
    x |= (x >> 1)
    x |= (x >> 2)
    # AND removes junk, we then widen x by multiplication and return
    # the inverse.
    x &= 0x1111111111111111
    x *= 15
    return ~x
end

"""
    create_nibble_mask(::Type{Gap}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent gaps.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Gap}, x::UInt64)
    return create_nibble_mask(x, 0x0000000000000000)
end

"""
    create_nibble_mask(::Type{Ambiguous}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent ambiguous sites.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Ambiguous}, x::UInt64)
    x = enumerate_nibbles(x)
    return create_nibble_mask(x, 0x2222222222222222) |
    create_nibble_mask(x, 0x3333333333333333) |
    create_nibble_mask(x, 0x4444444444444444)
end

"""
    create_nibble_mask(::Type{Pairdel}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites that should be
ignored, when counting pairwise mutations between sequences.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Pairdel}, x::UInt64)
    return create_nibble_mask(Gap, x) | create_nibble_mask(Ambiguous, x)
end

"""
    create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have a gap at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)
    return create_nibble_mask(Gap, a) | create_nibble_mask(Gap, b)
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
    create_nibble_mask(::Type{Pairdel}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have a nucleotide symbol that sould be ignored during
pairwise distance computation at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Pairdel}, a::UInt64, b::UInt64)
    return create_nibble_mask(Pairdel, a) | create_nibble_mask(Pairdel, b)
end

"""
    create_nibble_mask(::Type{Conserved}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides are the same.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Conserved}, a::UInt64, b::UInt64)
    pairdelmask = ~create_nibble_mask(Pairdel, a, b)
    a &= pairdelmask
    b &= pairdelmask
    return create_nibble_mask(a $ b, 0x0000000000000000) & pairdelmask
end

"""
    create_nibble_mask(::Type{Conserved}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides are different.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Mutated}, a::UInt64, b::UInt64)
    pairdelmask = ~create_nibble_mask(Pairdel, a, b)
    a &= pairdelmask
    b &= pairdelmask
    return ~create_nibble_mask(a $ b, 0x0000000000000000)
end


# Bitparallel site counting
# =========================

"""
    count_sites4(::Type{Gap}, x::UInt64)

Count the number of gap sites in a chunk of BioSequence{(DNA|RNA)Nucleotide{4}}
data.

Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Gap}, x::UInt64)
    return count_zero_nibbles(x)
end

"""
    count_sites4(::Type{Gap}, a::UInt64, b::UInt64)

Count the number of sites in two aligned chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data which contain gap characters.

Note that gap sites and empty unused segments of a UInt64 are both 0000, and so
furthur checking of this result would be required in higher level calling
functions.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Gap}, a::UInt64, b::UInt64)
    # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
    return count_sites4(Gap, a) + count_sites4(Gap, b) - count_sites4(Gap, a | b)
end

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
    count_sites4(::Type{Gap}, a::UInt64, b::UInt64)

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

"""
    count_sites4(::Type{Pairdel}, x::UInt64)

An _internal_ function _not for export_, which will count the number of sites in
a chunk of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored in
counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Pairdel}, x::UInt64)
    return count_zero_nibbles(~create_nibble_mask(Pairdel, x))
end

"""
    count_sites4(::Type{Pairdel}, x::UInt64)

An _internal_ function _not for export_, which will count the number of sites in
aligned chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data that would be ignored
in counts of mutations.
Such sites are defined as those with gaps or ambiguous characters in them.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Pairdel}, a::UInt64, b::UInt64)
    return count_zero_nibbles(~create_nibble_mask(Pairdel, a, b))
end

"""
    count_sites4(::Type{Conserved}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
Conserved between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
Conserved. For example, 'A' and 'R', or 'A' and '-' will not be counted.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Conserved}, a::UInt64, b::UInt64)
    pairdelmask = ~create_nibble_mask(Pairdel, a, b)
    a &= pairdelmask
    b &= pairdelmask
    diffs = a $ b
    initialZeros = count_zero_nibbles(pairdelmask)
    conservedZeros = count_zero_nibbles(diffs)
    return conservedZeros - initialZeros
end

"""
    count_sites4(::Type{Mutated}, a::UInt64, b::UInt64)

An _internal_ function, _not for export_, which will count the number of
any mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.

**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.

**This is an internal method and should not be exported.**
"""
@inline function count_sites4(::Type{Mutated}, a::UInt64, b::UInt64)
    #pairdelmask = ~create_nibble_mask(Pairdel, a, b)
    #a &= pairdelmask
    #b &= pairdelmask
    #diffs = a $ b
    #initialZeros = count_zero_nibbles(pairdelmask)
    #conservedZeros = count_zero_nibbles(diffs)
    #return 16 - (initialZeros + (conservedZeros - initialZeros))
    mutationMask = create_nibble_mask(Mutated, a, b)
    return count_one_nibbles(mutationMask)
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
    pairdelmask = ~create_nibble_mask(Pairdel, a, b)
    a &= pairdelmask
    b &= pairdelmask
    diffs = a $ b
    ts_filtered = diffs $ 0xAAAAAAAAAAAAAAAA
    return count_zero_nibbles(ts_filtered) + count_zero_nibbles(~ts_filtered)
end

@inline function count_sites4(::Type{Transversion}, a::UInt64, b::UInt64)
    pairdelmask = ~create_nibble_mask(Pairdel, a, b)
    a &= pairdelmask
    b &= pairdelmask
    diffs = a $ b
    tv_filtered = diffs $
    return count_zero_nibbles(tv_filtered) + count_zero_nibbles(~tv_filtered)
end
