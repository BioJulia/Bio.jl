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
    len = min(length(a), length(b))
    ritr = StepRange(width, step, len)
    width -= 1
    results = Vector{IntervalValue{out_type(T)}}(length(ritr))
    r = 1
    @inbounds for i in ritr
        idx = (i - width):i
        results[r] = IntervalValue(first(idx), last(idx), count(T, a[idx], b[idx]))
        r += 1
    end
    return results
end

function count_pairwise{T<:Site,N}(::Type{T}, seqs::Vararg{BioSequence,N})
    counts = Vector{out_type(T)}(Int((N * (N - 1)) / 2))
    c = 1
    @inbounds for i in 1:N, j in (i + 1):N
        counts[c] = count(T, seqs[i], seqs[j])
        c += 1
    end
    return PairwiseListMatrix(counts, false)
end
