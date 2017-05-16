# certain.jl
# ==========
#
# Define Certain sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

@inline function nibble_mask(::Type{Certain}, x::UInt64)
    return nibble_mask(enumerate_nibbles(x), 0x1111111111111111)
end

@inline function nibble_mask(::Type{Certain}, a::UInt64, b::UInt64)
    return nibble_mask(Certain, a) & nibble_mask(Certain, b)
end
