# issite methods for naieve site counting
# =======================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"Test whether a nucleotide site in a sequence is ambiguous."
issite{T<:Nucleotide}(::Type{Ambiguous}, a::T) = isambiguous(a)

"Test whether a nucleotide site of two aligned sequences has ambiguities."
@inline function issite{T<:Nucleotide}(::Type{Ambiguous}, a::T, b::T)
    return issite(Ambiguous, a) | issite(Ambiguous, b)
end

"Test whether a nucleotide site in a sequence is a gap character."
issite{T<:Nucleotide}(::Type{Gap}, a::T) = count_ones(a) == 0

"Test whether a nucleotide site of two aligned sequences has gap characters."
@inline function issite{T<:Nucleotide}(::Type{Gap}, a::T, b::T)
    return issite(Gap, a) | issite(Gap, b)
end

"Test whether a nucleotide site in a sequence is certain."
issite{T<:Nucleotide}(::Type{Certain}, a::T) = count_ones(a) == 1

"Test whether a nucleotide site of two aligned sequences has two certain nucleotides."
@inline function issite{T<:Nucleotide}(::Type{Certain}, a::T, b::T)
    return issite(Certain, a) & issite(Certain, b)
end

"Test whether a nucleotide site of two aligned sequences, constitutes a match."
issite{T<:Nucleotide}(::Type{Match}, a::T, b::T) = a == b

"Test whether a nucleotide site of two aligned sequences, constitutes a mismatch."
issite{T<:Nucleotide}(::Type{Mismatch}, a::T, b::T) = a != b

"Test whether a nucleotide site of two aligned seqences certainly constitutes a conserved site."
@inline function issite{T<:Nucleotide}(::Type{Conserved}, a::T, b::T)
    return issite(Certain, a, b) & issite(Match, a, b)
end

"Test whether a nucleotide site of two aligned seqences certainly constitutes a mutation."
@inline function issite{T<:Nucleotide}(::Type{Mutated}, a::T, b::T)
    return issite(Certain, a, b) & issite(Mismatch, a, b)
end

"Test whether a nucleotide site of two aligned sequences certainly constitutes a transition mutation."
@inline function issite{T<:Nucleotide}(::Type{Transition}, a::T, b::T)
    return issite(Mutated, a, b) & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

"Test whether a nucleotide site of two aligned sequences certainly constitutes a transversion mutation."
@inline function issite{T<:Nucleotide}(::Type{Transversion}, a::T, b::T)
    return issite(Mutated, a, b) & ((ispurine(a) & ispyrimidine(b)) | (ispyrimidine(a) & ispurine(b)))
end
