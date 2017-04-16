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

@inline function count_algorithm{A}(s::Mismatch, a::BioSequence{A}, b::BioSequence{A})
    return BITPAR
end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleotide site of two aligned sequences, constitutes a mismatch."
issite{T<:NucleicAcid}(::Type{Mismatch}, a::T, b::T) = a != b

# Methods for the bitparallel framework.
# --------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function count_bitpar(::Type{Mismatch}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_nonzero_nibbles(a $ b)
        end

        @inline function count_bitpar(::Type{Mismatch}, ::Type{$A{2}}, a::UInt64, b::UInt64)
            return count_nonzero_bitpairs(x $ y)
        end
    end
end
