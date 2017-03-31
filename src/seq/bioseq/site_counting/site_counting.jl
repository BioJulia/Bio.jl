# site_counting.jl
# ================
#
# A site counting framework for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

include("site_types.jl")
include("count_sites_naive.jl")

function Base.count{T<:Site}(::Type{T}, a::BioSequence)
    return count_sites_naieve(T, a)
end

function Base.count{T<:Site}(::Type{T}, a::BioSequence, b::BioSequence)
    return count_sites_naive(T, a, b)
end

function Base.count{T<:Site}(::Type{T}, a::BioSequence, b::BioSequence, width::Int, step::Int)
    ridx = 1:width
    len = min(length(a), length(b))
    ritr = StepRange(width, step, len)
    step -= 1
    results = Vector{IntervalValue{acc_type(T)}}(length(ritr))
    r = 1
    @inbounds for i in ritr
        idx = (i - step):i
        results[r] = IntervalValue(first(idx), last(idx), count(T, a[idx], b[idx]))
        r += 1
    end
    return results
end
