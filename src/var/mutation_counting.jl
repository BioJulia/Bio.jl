# phylo/mutation_counting.jl
# ==================
#
# Types and methods for counting different kinds of mutations between two
# DNA or RNA sequences.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


# Counting mutations
# ------------------

# Mutation types

abstract MutationType
immutable DifferentMutation <: MutationType end
immutable TransitionMutation <: MutationType end
immutable TransversionMutation <: MutationType end

"""
    is_ambiguous_strict{T<:Nucleotide}(a::T, b::T)

The strictest test of ambiguity, if either nucleotide is not
A, T, G, or C, then this function returns true.
"""
@inline function is_ambiguous_strict{T<:Nucleotide}(a::T, b::T)
    return isambiguous(a) | isambiguous(b)
end

"""
    is_mutation{T<:Nucleotide}(a::T, b::T, ::Type{DifferentMutation})

Test if two nucleotides constitute a `DifferentMutation`.
"""
@inline function is_mutation{T<:Nucleotide}(a::T, b::T, ::Type{DifferentMutation})
    return a != b
end

"""
    is_mutation{T<:Nucleotide}(a::T, b::T, ::Type{TransitionMutation})

Test if two nucleotides constitute a `TransitionMutation`.
"""
@inline function is_mutation{T<:Nucleotide}(a::T, b::T, ::Type{TransitionMutation})
    return (a != b) & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

"""
    is_mutation{T<:Nucleotide}(a::T, b::T, ::Type{TransversionMutation})

Test if two nucleotides constitute a `TransversionMutation`.
"""
@inline function is_mutation{T<:Nucleotide}(a::T, b::T, ::Type{TransversionMutation})
    return (a != b) & ((ispurine(a) & ispyrimidine(b)) | (ispyrimidine(a) & ispurine(b)))
end

"""
    count_mutations(sequences::Vector{BioSequence{A}}, ::Type{T})

Count the number of mutations between DNA sequences in a pairwise manner.

Different types of mutation can be counted:
`DifferentMutation`, `TransitionMutation`, `TransversionMutation`.

Returns a tuple of: 1. A vector containing the number of mutations between each,
possible pair of sequences, and 2. a vector containing the number of sites
considered (sites with any ambiguity characters are not considered) for each
possible pair of sequences.
"""
function count_mutations{A<:NucleotideAlphabets,T<:MutationType}(sequences::Vector{BioSequence{A}}, ::Type{T})
    # This method has been written with the aim of improving performance by taking
    # advantage of the memory layout of matrices of nucleotides, as well as
    # getting julia to emit simd code for the innermost loop.
    seqs = seqmatrix(sequences, :seq)
    S, N = size(seqs)
    c = binomial(N, 2)
    lengths = Vector{Int}(c)
    nmutations = Vector{Int}(c)
    target = 1
    @inbounds for i1 in 1:N
        for i2 in i1+1:N
            L = S
            Nd = 0
            for s in 1:S
                s1 = seqs[(i1 - 1) * S + s]
                s2 = seqs[(i2 - 1) * S + s]
                isamb = is_ambiguous_strict(s1, s2)
                ismut = is_mutation(s1, s2, T)
                L -= isamb
                Nd += !isamb & ismut
            end
            lengths[target] = L
            nmutations[target] = Nd
            target += 1
        end
    end
    return nmutations, lengths
end

"""
    count_mutations(sequences::Vector{BioSequence{A}}, ::Type{TransitionMutation}, ::Type{TransversionMutation})

Count the number of `TransitionMutation`s and `TransversionMutation`s in a
pairwise manner, between each possible pair of sequences.

Returns a tuple of: 1. A vector containing the number of transitions between
each, possible pair of sequences, and 2. a vector containing the number of
transversions between each, possible pair of sequences, and 3. a vector
containing the number of sites considered (sites with any ambiguity characters
are not considered) for each possible pair of sequences.
"""
function count_mutations{A<:NucleotideAlphabets}(sequences::Vector{BioSequence{A}}, ::Type{TransitionMutation}, ::Type{TransversionMutation})
    # This method has been written with the aim of improving performance by taking
    # advantage of the memory layout of matrices of nucleotides, as well as
    # getting julia to emit simd code for the innermost loop.
    seqs = seqmatrix(sequences, :seq)
    S, N = size(seqs)
    c = binomial(N, 2)
    lengths = Vector{Int}(c)
    ntransition = Vector{Int}(c)
    ntransversion = Vector{Int}(c)
    target = 1
    @inbounds for i1 in 1:N
        for i2 in i1+1:N
            L = S
            Ns = 0
            Nv = 0
            for s in 1:S
                s1 = seqs[(i1 - 1) * S + s]
                s2 = seqs[(i2 - 1) * S + s]
                isamb = is_ambiguous_strict(s1, s2)
                isdiff = s1 != s2
                istrans = (ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b))
                L -= isamb
                Ns += !isamb & isdiff & istrans
                Nv += !isamb & isdiff & !istrans
            end
            lengths[target] = L
            ntransition[target] = Ns
            ntransversion[target] = Nv
            target += 1
        end
    end
    return ntransition, ntransversion, lengths
end

function count_mutations{A<:NucleotideAlphabets}(sequences::Vector{BioSequence{A}}, ::Type{TransversionMutation}, ::Type{TransitionMutation})
    return count_mutations(sequences, TransitionMutation, TransversionMutation)
end
