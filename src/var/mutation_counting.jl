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

@inline function is_mutation{T<:Nucleotide}(a::T, b::T, t::Type{DifferentMutation})
    return a != b
end

@inline function is_mutation{T<:Nucleotide}(a::T, b::T, t::Type{TransitionMutation})
    return a != b & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

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
    count_mutations(seq1::BioSequence{A}, seq2::BioSequence{A}, t::Type{DifferentMutation})

Count the number of differences between two DNA sequences.

Returns a tuple of the number of sites that are different, and the number of sites
considered (sites with any ambiguity characters are not considered).
"""
function count_mutations{A<:NucleotideAlphabets}(seq1::BioSequence{A},
    seq2::BioSequence{A},
    t::Type{DifferentMutation})

    @assert length(seq1) == length(seq2)
    ndifferences = 0
    nsites = length(seq1)

    for (a, b) in zip(seq1, seq2)
        if is_ambiguous_strict(a, b)
            nsites -= 1
            continue
        end
        if a != b
            ndifferences += 1
        end
    end
    return ndifferences, nsites
end

"""
    count_mutations(seq1::BioSequence, seq2::BioSequence, t::Type{TransitionMutation})

Count the number of transition mutations between two DNA sequences.
Returns a tuple of the number of sites that are different, and the number of sites
considered (sites with any ambiguity characters are not considered).
"""
function count_mutations{A<:NucleotideAlphabets}(seq1::BioSequence{A},
    seq2::BioSequence{A},
    t::Type{TransitionMutation})

    @assert length(seq1) == length(seq2)
    ntransitions = 0
    nsites = length(seq1)

    for (a, b) in zip(seq1, seq2)
        if is_ambiguous_strict(a, b)
            nsites -= 1
            continue
        end
        if a != b && ((ispurine(a) && ispurine(b)) || (ispyrimidine(a) && ispyrimidine(b)))
            ntransitions += 1
        end
    end
    return ntransitions, nsites
end

"""
    count_mutations(seq1::BioSequence, seq2::BioSequence, t::Type{TransversionMutation})

Count the number of transversion mutations between two DNA sequences.
Returns a tuple of the number of sites that are different, and the number of sites
considered (sites with any ambiguity characters are not considered).
"""
function count_mutations{A<:NucleotideAlphabets}(seq1::BioSequence{A},
    seq2::BioSequence{A},
    t::Type{TransversionMutation})

    @assert length(seq1) == length(seq2)
    ntransversions = 0
    nsites = length(seq1)

    for (a, b) in zip(seq1, seq2)
        if is_ambiguous_strict(a, b)
            nsites -= 1
            continue
        end
        if a != b && ((ispurine(a) && ispyrimidine(b)) || (ispyrimidine(a) && ispurine(b)))
            ntransversions += 1
        end
    end
    return ntransversions, nsites
end

"""
    count_mutations(seq1::BioSequence, seq2::BioSequence, t1::Type{TransitionMutation}, t2::Type{TransversionMutation})

Count the number of transition and transversion mutations between two DNA sequences.
"""
function count_mutations{A<:NucleotideAlphabets}(seq1::BioSequence{A},
    seq2::BioSequence{A},
    t1::Type{TransitionMutation},
    t2::Type{TransversionMutation})

    @assert length(seq1) == length(seq2)
    ndifferences = 0
    ntransitions = 0
    nsites = length(seq1)

    for (a, b) in zip(seq1, seq2)
        if is_ambiguous_strict(a, b)
            nsites -= 1
            continue
        end
        if a != b
            ndifferences += 1
            if (ispurine(a) && ispurine(b)) || (ispyrimidine(a) && ispyrimidine(b))
                ntransitions += 1
            end
        end
    end
    ntransversions = ndifferences - ntransitions
    return ntransitions, ntransversions, nsites
end

function count_mutations{A<:NucleotideAlphabets}(a::BioSequence{A},
    b::BioSequence{A},
    t1::Type{TransversionMutation},
    t2::Type{TransitionMutation})

    return count_mutations(a, b, t2, t1)
end
