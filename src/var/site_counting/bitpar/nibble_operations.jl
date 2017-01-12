# Nibble operations
# =================
#
# Bit operations for working on nibbles - bases in 4 bit encoded nucleotides.
#
# These functions probably should not be exported, as it's easy to do things
# very wrong with them if you are not careful. However, they allow the
# identification of sites to be fast, by taking advantage of the
# BioSequence compressed representation and avoiding memory expensive
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
(In the BioSequence 4 bit encoding, each nibble is one RNA or DNA nucleotide).

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
    # XOR with the desired values. So matching nibbles will be 0000.
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


## Indel sites

"""
    create_nibble_mask(::Type{Indel}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent gaps.

Care should be taken as unused nibbles, as well as indels will be allowed through
this mask.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Indel}, x::UInt64)
    return create_nibble_mask(x, 0x0000000000000000)
end

"""
    create_nibble_mask(::Type{Indel}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have a gap at a given nibble.

Care should be taken as unused nibbles, as well as gaps will be allowed through
this mask.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Indel}, a::UInt64, b::UInt64)
    return create_nibble_mask(Indel, a) | create_nibble_mask(Indel, b)
end


## Certain sites

"""
    create_nibble_mask(::Type{Certain}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites that should be
considered, when counting pairwise mutations between sequences.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Certain}, x::UInt64)
    return create_nibble_mask(enumerate_nibbles(x), 0x1111111111111111)
end

"""
    create_nibble_mask(::Type{Certain}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites that should be
considered during pairwise distance computation at a given nibble.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Certain}, a::UInt64, b::UInt64)
    return create_nibble_mask(Certain, a) & create_nibble_mask(Certain, b)
end


## Ambiguous sites

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


## Matching sites

"""
    create_nibble_mask(::Type{Match}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the biological
symbols are the same.

Care should be taken as unused sites, as well as indels may be allowed through
this mask as both are seen as a pair of matching 0000 nibbles.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Match}, a::UInt64, b::UInt64)
    return create_nibble_mask(0x0000000000000000, a $ b)
end


## Mismatching sites

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


## Conserved sites

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


## Mutated sites

"""
    create_nibble_mask(::Type{Mutated}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides are different.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Mutated}, a::UInt64, b::UInt64)
    certainmask = create_nibble_mask(Certain, a, b)
    return (~create_nibble_mask(a $ b, 0x0000000000000000)) & certainmask
end

"""
    create_nibble_mask(::Type{Transition}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides constitute a transition mutation.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Transition}, a::UInt64, b::UInt64)
    ts_filter = (a $ b) & create_nibble_mask(Certain, a, b)
    return create_nibble_mask(ts_filter, 0xAAAAAAAAAAAAAAAA) |
           create_nibble_mask(ts_filter, 0x5555555555555555)
end

"""
    create_nibble_mask(::Type{Transversion}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides constitute a transition mutation.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Transversion}, a::UInt64, b::UInt64)
    diffs = (a $ b) & create_nibble_mask(Certain, a, b)
    return create_nibble_mask(diffs, 0x3333333333333333) |
           create_nibble_mask(diffs, 0x6666666666666666) |
           create_nibble_mask(diffs, 0xCCCCCCCCCCCCCCCC) |
           create_nibble_mask(diffs, 0x9999999999999999)

end
