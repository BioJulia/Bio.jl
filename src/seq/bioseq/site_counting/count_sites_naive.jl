# count_sites_naive
# =================
#
# Methods for site counting in sequences using a naive algorithm.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

NC = NaiveCount

@inline function Base.count{S<:Site,A<:Alphabet}(::Type{S}, ::Type{NC}, a::BioSequence{A})
    k = start_counter(S, A)
    @inbounds for x in eachindex(a)
        k = update_counter(k, issite(S, a, x))
    end
    return k
end

@inline function Base.count{S<:Site,A<:Alphabet,B<:Alphabet}(::Type{S}, ::Type{NC}, a::BioSequence{A}, b::BioSequence{B})
    k = start_counter(S, A, B)
    @inbounds for idx in 1:min(endof(a), endof(b))
        k = update_counter(k, issite(S, a, b, idx))
    end
    return k
end
