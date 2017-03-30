# count_sites_naive
# =================
#
# Methods for site counting in sequences using a naive method.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Methods for a single sequence.

function count_sites_naive{T<:Site}(::Type{T}, a::BioSequence)
    k = start_count(T)
    @inbounds for x in a
         k = update_count(k..., issite(T, x)...)
    end
    return k
end

function count_sites_naive{T<:Site}(::Type{T}, a::BioSequence, b::BioSequence)
    k = start_count(T)
    @inbounds for idx in 1:min(endof(a), endof(b))
        k = update_count(k..., issite(T, a[idx], b[idx])...)
    end
    return k
end
