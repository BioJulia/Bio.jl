# Types and methods for counting mutations
# ========================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

#=typealias FourBitAlphs Union{DNAAlphabet{4},RNAAlphabet{4}}
typealias TwoBitAlphs Union{DNAAlphabet{2},RNAAlphabet{2}}
typealias NucAlphs Union{DNAAlphabet, RNAAlphabet}

"""
Count the number of mutations between nucleotide sequences in a pairwise manner,
using a sliding window of a given width and step.
"""
function count_sites{T<:Site,A<:Alphabet}(::Type{T}, sequences::Array{BioSequence{A}}, width::Int, step::Int)
    len = length(sequences)
    counts = PairwiseListMatrix(Vector{IntervalValue{Int}}, len, false)
    @inbounds for i in 1:len, j in (i + 1):len
        counts[i, j] = count_sites(T, sequences[i], sequences[j], width, step)
    end
    return counts
end

"""
Count the number of mutations between nucleotide sequences in a pairwise manner.
"""
function count_sites{T<:Mutation,A<:Alphabet}(::Type{T}, sequences::Array{BioSequence{A}})
    len = length(sequences)
    counts = PairwiseListMatrix(Int, len, false)
    undetermined = PairwiseListMatrix(Int, len, false)
    @inbounds for i in 1:len, j in (i + 1):len
        counts[i, j], undetermined[i, j] = count_sites(T, sequences[i], sequences[j])
    end
    return counts, undetermined
end

"""
Count the number of mutations between nucleotide sequences in a pairwise manner,
using a sliding window of a given width and step.
"""
function count_sites{T<:Mutation,A<:Alphabet}(::Type{T}, sequences::Array{BioSequence{A}}, width::Int, step::Int)
    len = length(sequences)
    counts = PairwiseListMatrix(Vector{IntervalValue{Int}}, len, false)
    undetermined = PairwiseListMatrix(Vector{IntervalValue{Int}}, len, false)
    @inbounds for i in 1:len, j in (i + 1):len
        counts[i, j], undetermined[i, j] = count_sites(T, sequences[i], sequences[j], width, step)
    end
    return counts, undetermined
end

"Count the number of sites of a given SiteType `T` in a biological sequence."
function count_sites{T<:Site}(::Type{T}, seq::BioSequence)
    return count_sites_naive(T, seq)
end
=#
