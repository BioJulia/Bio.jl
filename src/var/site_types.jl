# Extending BioSequence site counting to include mutations
# ========================================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

Base.zero(::Type{Tuple{Int,Int}}) = zero(Int), zero(Int)


abstract Mutation <: Site
start_count{T<:Mutation}(::Type{T}) = 0, 0
out_type{T<:Mutation}(::Type{T}) = Tuple{Int,Int}
update_count(kacc::Int, dacc::Int, k::Bool, d::Bool) = kacc + k, dacc + d

@inline function issite{T<:NucleicAcid,M<:Mutation}(::Type{M}, a::T, b::T)
    k = ischange(M, a, b)
    d = issite(Certain, a, b)
    return k & d, !d
end

"""
A `Mismatch` site describes a site where two aligned nucleotides are definately
conserved. By definately conserved this means that the symbols of the site are
non-ambiguity symbols, and they are the same symbol.
"""
immutable Conserved <: Mutation end

@inline function ischange{T<:NucleicAcid}(::Type{Conserved}, a::T, b::T)
    return issite(Match, a, b)
end

"""
A `Mutated` site describes a site where two aligned nucleotides are definately
mutated. By definately mutated this means that the symbols of the site are
non-ambiguity symbols, and they are not the same symbol.
"""
immutable Mutated <: Mutation end

@inline function ischange{T<:NucleicAcid}(::Type{Mutated}, a::T, b::T)
    return issite(Mismatch, a, b)
end

"""
A `Transition` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transition mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transition mutation: i.e. A<->G, or C<->T.
"""
immutable Transition <: Mutation end

@inline function ischange{T<:NucleicAcid}(::Type{Transition}, a::T, b::T)
    return issite(Mismatch, a, b) & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

"""
A `Transversion` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transversion mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transversion mutation: i.e. A<->C,
A<->T, G<->T, G<->C.
"""
immutable Transversion <: Mutation end

@inline function ischange{T<:NucleicAcid}(::Type{Transversion}, a::T, b::T)
    return issite(Mismatch, a, b) & ((ispurine(a) & ispyrimidine(b)) | (ispyrimidine(a) & ispurine(b)))
end
