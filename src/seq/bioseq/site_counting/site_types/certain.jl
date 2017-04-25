# certain.jl
# ==========
#
# Define Certain sites for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol or a gap.
"""
immutable Certain <: Site end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleic acid in a sequence is certain."
issite(::Type{Certain}, a::BioSequence, idx) = iscertain(a[idx])

"Test whether a nucleotide site of two aligned sequences has two certain NucleicAcids."
@inline function issite(::Type{Certain}, a::BioSequence, b::BioSequence, idx)
    return issite(Certain, a, idx) & issite(Certain, b, idx)
end

# Methods for the bitparallel framework.
# --------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function count_algorithm(s::Certain, a::BioSequence{$A{2}}, b::BioSequence{$A{2}})
            return ALL
        end

        @inline function count_algorithm(s::Certain, a::BioSequence{$A{4}}, b::BioSequence{$A{4}})
            return BITPAR
        end

        @inline function count_bitpar(::Type{Certain}, ::Type{$A{4}}, x::UInt64)
            x = enumerate_nibbles(x)
            x $= 0x1111111111111111
            return count_zero_nibbles(x)
        end

        @inline function count_bitpar(::Type{Certain}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            x = enumerate_nibbles(a) $ 0x1111111111111111
            y = enumerate_nibbles(b) $ 0x1111111111111111
            return count_zero_nibbles(x | y)
        end
    end
end

@inline function nibble_mask(::Type{Certain}, x::UInt64)
    return nibble_mask(enumerate_nibbles(x), 0x1111111111111111)
end

@inline function nibble_mask(::Type{Certain}, a::UInt64, b::UInt64)
    return nibble_mask(Certain, a) & nibble_mask(Certain, b)
end
