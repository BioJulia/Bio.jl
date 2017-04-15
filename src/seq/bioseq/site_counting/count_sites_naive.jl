# count_sites_naive
# =================
#
# Methods for site counting in sequences using a naive algorithm.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

NC = NaiveCount

@inline function Base.count{T<:Site}(site::T, alg::NC, a::BioSequence)
    k = start_counter(site, a)
    @inbounds for x in a
        k = update_counter(k..., issite(T, x)...)
    end
    return k
end

@inline function Base.count{T<:Site}(site::T, alg::NC, a::BioSequence, b::BioSequence)
    k = start_counter(site, a, b)
    @inbounds for idx in 1:min(endof(a), endof(b))
        k = update_counter(k..., issite(T, a[idx], b[idx])...)
    end
    return k
end
