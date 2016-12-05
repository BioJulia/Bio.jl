# phylo/mutation_counting.jl
# ==========================
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

"""
`AnyMutation` describes a site where two aligned nucleotides are not the
same.

Every kind of difference is counted.
"""
immutable AnyMutation <: MutationType end

"""
`TransitionMutation` describes a situation with two aligned nucleotides, where a
purine has mutated into another purine, or a pyrimadine has mutated into another
pyrimadine.

Possible transition mutations are:
A <-> G
C <-> T
"""
immutable TransitionMutation <: MutationType end

"""
`TransversionMutation` describes a situation with two aligned nucleotides,
where a purine has mutated into a pyrimadine or vice versa.

Possible transversion mutations are:
A <-> T
C <-> G
A <-> C
T <-> G
"""
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
    is_mutation{T<:Nucleotide}(::Type{AnyMutation}, a::T, b::T)

Test if two nucleotides constitute a `AnyMutation`.
"""
@inline function is_mutation{T<:Nucleotide}(::Type{AnyMutation}, a::T, b::T)
    return a != b
end

"""
    is_mutation{T<:Nucleotide}(::Type{TransitionMutation}, a::T, b::T)

Test if two nucleotides constitute a `TransitionMutation`.
"""
@inline function is_mutation{T<:Nucleotide}(::Type{TransitionMutation}, a::T, b::T)
    return (a != b) & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

"""
    is_mutation{T<:Nucleotide}(::Type{TransversionMutation}, a::T, b::T)

Test if two nucleotides constitute a `TransversionMutation`.
"""
@inline function is_mutation{T<:Nucleotide}(::Type{TransversionMutation}, a::T, b::T)
    return (a != b) & ((ispurine(a) & ispyrimidine(b)) | (ispyrimidine(a) & ispurine(b)))
end





"""
    flagmutations{M<:MutationType,N<:Nucleotide}(::Type{M}, seqs::Matrix{N})

For every pair of sequences, flag which base positions are mutations of type `T`.

This function returns a tuple of two matrices. Both matrices contain logical
values, and contain a column for every pairwise combination of sequences,
and a number of rows equal to the length of the aligned sequences.

For the first matrix, for a cell of a given column,
a value of `true` indicates that there is a mutation of type `T` between two
sequences at that base position.

For the second matrix, for a cell of a given column,
a value of `true` indicates that there is an ambiguous base in either of the two
sequences at that base position, consequently it is not possible to conclude
whether there is a mutation at that base position or not.

**Note: This method assumes that the sequences are stored in the `Matrix{N}`
provided as `seqs` in sequence major order i.e. each column of the matrix is one
complete nucleotide sequence.**

# Example

```julia
dnas = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
2-element Array{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},1}:
 ATTG-ACCTGGNTTTCCGAA
 A-ACAGAGTATACRGTCGTC

julia> m = seqmatrix(dnas, :seq)
20×2 Array{Bio.Seq.DNANucleotide,2}:
  DNA_A    DNA_A
  DNA_T    DNA_Gap
  DNA_T    DNA_A
  DNA_G    DNA_C
  DNA_Gap  DNA_A
  DNA_A    DNA_G
  DNA_C    DNA_A
  DNA_C    DNA_G
  DNA_T    DNA_T
  DNA_G    DNA_A
  DNA_G    DNA_T
  DNA_N    DNA_A
  DNA_T    DNA_C
  DNA_T    DNA_R
  DNA_T    DNA_G
  DNA_C    DNA_T
  DNA_C    DNA_C
  DNA_G    DNA_G
  DNA_A    DNA_T
  DNA_A    DNA_C

julia> r = flagmutations(AnyMutation, m)
(
Bool[false; false; … ; true; true],

Bool[false; true; … ; false; false])

julia> r[1]
20×1 Array{Bool,2}:
 false
 false
  true
  true
 false
  true
  true
  true
 false
  true
  true
 false
  true
 false
  true
  true
 false
 false
  true
  true
```

In the above example of two sequences, positions 3:4, 7:8 are mutations and
positions 1:2, 5, 19:20 are not mutations.

"""
function flagmutations{M<:MutationType,N<:Nucleotide}(::Type{M}, seqs::Matrix{N})
    seqsize, nseqs = size(seqs)
    ismutant = Matrix{Bool}(seqsize, binomial(nseqs, 2))
    isambiguous = Matrix{Bool}(seqsize, binomial(nseqs, 2))
    col = 1
    @inbounds for i1 in 1:nseqs
        for i2 in i1+1:nseqs
            for s in 1:seqsize
                s1 = seqs[(i1 - 1) * seqsize + s]
                s2 = seqs[(i2 - 1) * seqsize + s]
                isamb = is_ambiguous_strict(s1, s2)
                ismut = is_mutation(M, s1, s2)
                isambiguous[(col - 1) * seqsize + s] = isamb
                ismutant[(col - 1) * seqsize + s] = !isamb & ismut
            end
            col += 1
        end
    end
    return ismutant, isambiguous
end

function flagmutations{M<:MutationType,A<:NucleotideAlphabet}(::Type{M}, seqs::Vector{BioSequence{A}})
    return flagmutations(M, seqmatrix(seqs, :seq))
end

"""
    count_mutations{T<:MutationType,N<:Nucleotide}(::Type{T}, seqs::Matrix{N})

Count the number of mutations between DNA sequences in a pairwise manner.

Different types of mutation can be counted:
`AnyMutation`, `TransitionMutation`, `TransversionMutation`.

Returns a tuple of: 1. A vector containing the number of mutations between each,
possible pair of sequences, and 2. a vector containing the number of sites
considered (sites with any ambiguity characters are not considered) for each
possible pair of sequences.

**Note: This method assumes that the sequences are stored in the `Matrix{N}`
provided as `seqs` in sequence major order i.e. each column of the matrix is one
complete nucleotide sequence.**
"""
function count_mutations{T<:MutationType,N<:Nucleotide}(::Type{T}, seqs::Matrix{N})
    # This method has been written with the aim of improving performance by taking
    # advantage of the memory layout of matrices of nucleotides, as well as
    # getting julia to emit simd code for the innermost loop.
    seqsize, nseqs = size(seqs)
    c = binomial(nseqs, 2)
    lengths = Vector{Int}(c)
    nmutations = Vector{Int}(c)
    target = 1
    @inbounds for i1 in 1:nseqs
        for i2 in i1+1:nseqs
            L = seqsize
            Nd = 0
            for s in 1:seqsize
                s1 = seqs[(i1 - 1) * seqsize + s]
                s2 = seqs[(i2 - 1) * seqsize + s]
                isamb = is_ambiguous_strict(s1, s2)
                ismut = is_mutation(T, s1, s2)
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
    count_mutations{T<:MutationType,A<:NucleotideAlphabet}(::Type{T}, sequences::Vector{BioSequence{A}})

Count the number of mutations between DNA sequences in a pairwise manner.

Different types of mutation can be counted:
`AnyMutation`, `TransitionMutation`, `TransversionMutation`.

Returns a tuple of: 1. A vector containing the number of mutations between each,
possible pair of sequences, and 2. a vector containing the number of sites
considered (sites with any ambiguity characters are not considered) for each
possible pair of sequences.
"""
function count_mutations{T<:MutationType,A<:NucleotideAlphabet}(::Type{T}, sequences::Vector{BioSequence{A}})
    seqs = seqmatrix(sequences, :seq)
    return count_mutations(T, seqs)
end


function flagmutations{N<:Nucleotide}(::Type{TransitionMutation}, ::Type{TransversionMutation}, seqs::Matrix{N})
    seqsize, nseqs = size(seqs)
    istransition = Matrix{Bool}(seqsize, binomial(nseqs, 2))
    istransversion = Matrix{Bool}(seqsize, binomial(nseqs, 2))
    isambiguous = Matrix{Bool}(seqsize, binomial(nseqs, 2))
    col = 1
    @inbounds for i1 in 1:nseqs
        s1offset = (i1 - 1) * seqsize
        for i2 in i1+1:nseqs
            s2offset = (i2 - 1) * seqsize
            resoffset = (col - 1) * seqsize
            for s in 1:seqsize
                s1 = seqs[s1offset + s]
                s2 = seqs[s2offset + s]
                isamb = is_ambiguous_strict(s1, s2)
                isdiff = s1 != s2
                ists = is_mutation(TransitionMutation, s1, s2)
                isambiguous[resoffset + s] = isamb
                istransition[resoffset + s] = !isamb & isdiff & ists
                istransversion[resoffset + s] = !isamb & isdiff & !ists
            end
            col += 1
        end
    end
    return istransition, istransversion, isambiguous
end

function flagmutations{A<:NucleotideAlphabet}(::Type{TransitionMutation}, ::Type{TransversionMutation}, seqs::Vector{BioSequence{A}})
    return flagmutations(TransitionMutation, TransversionMutation, seqmatrix(seqs, :seq))
end



"""
    count_mutations{N<:Nucleotide}(::Type{TransitionMutation}, ::Type{TransversionMutation}, sequences::Matrix{N})

Count the number of `TransitionMutation`s and `TransversionMutation`s in a
pairwise manner, between each possible pair of sequences.

Returns a tuple of: 1. A vector containing the number of transitions between
each, possible pair of sequences, and 2. a vector containing the number of
transversions between each, possible pair of sequences, and 3. a vector
containing the number of sites considered (sites with any ambiguity characters
are not considered) for each possible pair of sequences.

**Note: This method assumes that the sequences are stored in the `Matrix{N}`
provided as `seqs` in sequence major order i.e. each column of the matrix is one
complete nucleotide sequence.**
"""
function count_mutations{A<:Nucleotide}(::Type{TransitionMutation}, ::Type{TransversionMutation}, seqs::Matrix{A})
    # This method has been written with the aim of improving performance by taking
    # advantage of the memory layout of matrices of nucleotides, as well as
    # getting julia to emit simd code for the innermost loop.
    seqsize, nseqs = size(seqs)
    c = binomial(nseqs, 2)
    lengths = Vector{Int}(c)
    ntransition = Vector{Int}(c)
    ntransversion = Vector{Int}(c)
    target = 1
    @inbounds for i1 in 1:nseqs
        for i2 in i1+1:nseqs
            L = seqsize
            Ns = 0
            Nv = 0
            for s in 1:seqsize
                s1 = seqs[(i1 - 1) * seqsize + s]
                s2 = seqs[(i2 - 1) * seqsize + s]
                isamb = is_ambiguous_strict(s1, s2)
                isdiff = s1 != s2
                istrans = (ispurine(s1) & ispurine(s2)) | (ispyrimidine(s1) & ispyrimidine(s2))
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

function count_mutations{A<:Nucleotide}(::Type{TransversionMutation}, ::Type{TransitionMutation}, seqs::Matrix{A})
    return count_mutations(TransitionMutation, TransversionMutation, seqs)
end

"""
    count_mutations{A<:NucleotideAlphabet}(::Type{TransitionMutation}, ::Type{TransversionMutation}, sequences::Vector{BioSequence{A}})

Count the number of `TransitionMutation`s and `TransversionMutation`s in a
pairwise manner, between each possible pair of sequences.

Returns a tuple of: 1. A vector containing the number of transitions between
each, possible pair of sequences, and 2. a vector containing the number of
transversions between each, possible pair of sequences, and 3. a vector
containing the number of sites considered (sites with any ambiguity characters
are not considered) for each possible pair of sequences.
"""
function count_mutations{A<:NucleotideAlphabet}(::Type{TransitionMutation}, ::Type{TransversionMutation}, sequences::Vector{BioSequence{A}})
    seqs = seqmatrix(sequences, :seq)
    return count_mutations(TransitionMutation, TransversionMutation, seqs)
end

function count_mutations{A<:NucleotideAlphabet}(::Type{TransversionMutation}, ::Type{TransitionMutation}, sequences::Vector{BioSequence{A}})
    return count_mutations(TransitionMutation, TransversionMutation, sequences)
end
