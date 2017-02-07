# issite methods for naieve site counting
# =======================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"Test whether a nucleic acid in a sequence is ambiguous."
issite{T<:NucleicAcid}(::Type{Ambiguous}, a::T) = isambiguous(a)

"Test whether a nucleotide site of two aligned sequences has ambiguities."
@inline function issite{T<:NucleicAcid}(::Type{Ambiguous}, a::T, b::T)
    return issite(Ambiguous, a) | issite(Ambiguous, b)
end

"Test whether a nucleic acid in a sequence is a gap character."
issite{T<:NucleicAcid}(::Type{Gap}, a::T) = isgap(a)

"Test whether a nucleotide site of two aligned sequences has gap characters."
@inline function issite{T<:NucleicAcid}(::Type{Gap}, a::T, b::T)
    return issite(Gap, a) | issite(Gap, b)
end

"Test whether a nucleic acid in a sequence is certain."
issite{T<:NucleicAcid}(::Type{Certain}, a::T) = iscertain(a)

"Test whether a nucleotide site of two aligned sequences has two certain NucleicAcids."
@inline function issite{T<:NucleicAcid}(::Type{Certain}, a::T, b::T)
    return issite(Certain, a) & issite(Certain, b)
end

"Test whether a nucleotide site of two aligned sequences, constitutes a match."
issite{T<:NucleicAcid}(::Type{Match}, a::T, b::T) = a == b

"Test whether a nucleotide site of two aligned sequences, constitutes a mismatch."
issite{T<:NucleicAcid}(::Type{Mismatch}, a::T, b::T) = a != b

"Test whether a nucleotide site of two aligned seqences certainly constitutes a conserved site."
@inline function issite{T<:NucleicAcid}(::Type{Conserved}, a::T, b::T)
    return issite(Certain, a, b) & issite(Match, a, b)
end

"Test whether a nucleotide site of two aligned seqences certainly constitutes a mutation."
@inline function issite{T<:NucleicAcid}(::Type{Mutated}, a::T, b::T)
    return issite(Certain, a, b) & issite(Mismatch, a, b)
end

"Test whether a nucleotide site of two aligned sequences certainly constitutes a transition mutation."
@inline function issite{T<:NucleicAcid}(::Type{Transition}, a::T, b::T)
    return issite(Mutated, a, b) & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

"Test whether a nucleotide site of two aligned sequences certainly constitutes a transversion mutation."
@inline function issite{T<:NucleicAcid}(::Type{Transversion}, a::T, b::T)
    return issite(Mutated, a, b) & ((ispurine(a) & ispyrimidine(b)) | (ispyrimidine(a) & ispurine(b)))
end
