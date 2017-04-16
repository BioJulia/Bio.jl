# site_counting.jl
# ================
#
# A site counting framework for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract CountAlgorithm
immutable NaiveCount <: CountAlgorithm end
immutable BitparCount <: CountAlgorithm end
immutable NoneCount <: CountAlgorithm end
immutable AllCount <: CountAlgorithm end

include("binary_operations.jl")
include("site_types/site_types.jl")
include("count_sites_naive.jl")
include("count_sites_bitpar.jl")

@inline function Base.count{S<:Site}(::Type{S}, a::BioSequence)
    return count(S, NaiveCount, a)
end

@inline function Base.count{S<:Site}(::Type{S}, a::BioSequence, b::BioSequence)
    return count(S, NaiveCount, a, b)
end

function Base.count{S<:Site,A<:Alphabet,B<:Alphabet}(::Type{S}, a::BioSequence{A}, b::BioSequence{B}, width::Int, step::Int)
    len = min(length(a), length(b))
    ritr = StepRange(width, step, len)
    width -= 1
    results = Vector{IntervalValue{Int,counter_type(S, A, B)}}(length(ritr))
    r = 1
    @inbounds for i in ritr
        idx = (i - width):i
        results[r] = IntervalValue(first(idx), last(idx), count(S, a[idx], b[idx]))
        r += 1
    end
    return results
end

function count_pairwise{S<:Site,N}(::Type{S}, seqs::Vararg{BioSequence,N})
    counts = Vector{counter_type(S, a, b)}(Int((N * (N - 1)) / 2))
    c = 1
    @inbounds for i in 1:N, j in (i + 1):N
        counts[c] = count(S, seqs[i], seqs[j])
        c += 1
    end
    return PairwiseListMatrix(counts, false)
end
