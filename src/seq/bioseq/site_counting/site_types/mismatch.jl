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

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        count_algorithm(::Type{Mismatch}, a::BioSequence{$A{4}}, b::BioSequence{$A{4}}) = BITPAR
        count_algorithm(::Type{Mismatch}, a::BioSequence{$A{2}}, b::BioSequence{$A{2}}) = BITPAR
    end
end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleotide site of two aligned sequences, constitutes a mismatch."
issite(::Type{Mismatch}, a::BioSequence, b::BioSequence, idx) = a[idx] != b[idx]

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
