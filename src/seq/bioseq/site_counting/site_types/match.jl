# match.jl
# ========
#
# Define Matches for the site-counting framework.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
immutable Match <: Site end

@inline function count_algorithm{A<:Union{DNAAlphabet,RNAAlphabet}}(s::Match, a::BioSequence{A}, b::BioSequence{A})
    return BITPAR
end

# Methods for the naive framework.
# --------------------------------

"Test whether a nucleotide site of two aligned sequences, constitutes a match."
issite{T<:NucleicAcid}(::Type{Match}, a::T, b::T) = a == b

# Methods for the bitparallel framework.
# --------------------------------------

@inline correct_endspace{A<:Alphabet}(::Type{Match}, ::Type{A}) = true

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        @inline function count_bitpar(::Type{Match}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_zero_nibbles(a $ b)
        end

        @inline function count_bitpar(::Type{Match}, ::Type{$A{2}}, a::UInt64, b::UInt64)
            return count_zero_bitpairs(a $ b)
        end
    end
end
