# gap.jl
# =============
#
# Define Gaps for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
An `Gap` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
immutable Gap <: Site end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleic acid in a sequence is a gap character."
issite{T<:NucleicAcid}(::Type{Gap}, a::T) = isgap(a)

"Test whether a nucleotide site of two aligned sequences has gap characters."
@inline function issite{T<:NucleicAcid}(::Type{Gap}, a::T, b::T)
    return issite(Gap, a) | issite(Gap, b)
end

# Methods for the bitparallel framework.
# --------------------------------------

@inline correct_endspace(::Type{Gap}) = true

"""
    create_nibble_mask(::Type{Gap}, x::UInt64)

Create a mask of the nibbles in a chunk of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent gaps.

Care should be taken as unused nibbles, as well as indels will be allowed through
this mask.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Gap}, x::UInt64)
    return create_nibble_mask(x, 0x0000000000000000)
end

"""
    create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where either or
both of the two chunks have a gap at a given nibble.

Care should be taken as unused nibbles, as well as gaps will be allowed through
this mask.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Gap}, a::UInt64, b::UInt64)
    return create_nibble_mask(Gap, a) | create_nibble_mask(Gap, b)
end

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
