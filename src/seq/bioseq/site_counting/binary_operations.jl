# Binary operations
# =================
#
# Bit operations for working on nibbles (bases in 4 bit encoded format)
# and bit-pairs (bases in 2 bit encoded format).
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

@inline function count_nonzero_nibbles(x::UInt64)
    out = UInt64(0)
    out |= x & 0x1111111111111111
    out |= (x & 0x2222222222222222) >> 1
    out |= (x & 0x4444444444444444) >> 2
    out |= (x & 0x8888888888888888) >> 3
    return count_ones(out)
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
    return 16 - count_nonzero_nibbles(x)
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
    out = x & 0x1111111111111111
    out &= (x & 0x2222222222222222) >> 1
    out &= (x & 0x4444444444444444) >> 2
    out &= (x & 0x8888888888888888) >> 3
    return count_ones(out)
end

# Nibble masking functions
# ------------------------

"""
    nibble_mask(x::UInt64, value::UInt64)

Create a mask for the nibbles (groups of four bits) in a 64 bit integer `x`
that match a given value dictated by the pattern in `value`.

**This is an internal method and should not be exported.**
"""
@inline function nibble_mask(value::UInt64, x::UInt64)
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

# Bit-pair operations
# ===================

@inline function count_nonzero_bitpairs(x::UInt64)
    out = x & 0x5555555555555555
    out |= (x & 0xAAAAAAAAAAAAAAAA) >> 1
    return count_ones(out)
end

@inline function count_zero_bitpairs(x::UInt64)
    return 32 - count_nonzero_bitpairs(x)
end

@inline function count_one_bitpairs(x::UInt64)
    out = x & 0x5555555555555555
    out &= (x & 0xAAAAAAAAAAAAAAAA) >> 1
    return count_ones(out)
end
