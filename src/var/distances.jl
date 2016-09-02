# var/distances.jl
# ==================
#
# Types and methods for computing evolutionary and genetic distances.
#
# Part of the Bio.Var module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Types
# -----

abstract EvolutionaryDistance
abstract UncorrectedDistance <: EvolutionaryDistance
abstract CorrectedDistance <: EvolutionaryDistance
abstract TsTv <: CorrectedDistance

"""
A distance which is the count of the mutations of type T that exist between the
two sequences.
"""
immutable Count{T} <: UncorrectedDistance end

"""
Proportion{T} is a distance which is the count of the mutations of type T that
exist between the two biological sequences, divided by the number of valid sites
examined (sites which don't have gap or ambiguous symbols).

In other words this so called p-distance is simply the proportion of sites
between each pair of sequences, that are mutated (again where T determines
what kind of mutation).
"""
immutable Proportion{T} <: UncorrectedDistance end

"""
The JukesCantor69 distance is a p-distance adjusted/corrected by the
substitution model developed by Jukes and Cantor in 1969.

The Jukes and Cantor model assumes that all substitutions
(i.e. a change of a base by another one) have the same probability.
This probability is the same for all sites along the DNA sequence.
"""
immutable JukesCantor69 <: CorrectedDistance end

"""
The Kimura80 distance uses a substitution model developed by Kimura in 1980.
It is somtimes called Kimura's 2 parameter distance.

The model makes the same assumptions as Jukes and Cantor's model, but with a
crucial difference: two-kinds of mutation are considered
called Transitions and Transversions.
Transitions and transversions can occur with different probabilities in this
model, however, both transition and transversion rates/probabilities are the
same for all sites along the DNA sequence.
"""
immutable Kimura80 <: TsTv end


# Distance computation internals
# ------------------------------

## Jukes and Cantor 1969 distance computation.

@inline function expected_distance(::Type{JukesCantor69}, p::Float64)
    return -0.75 * log(1 - 4 * p / 3)
end

@inline function variance(::Type{JukesCantor69}, p::Float64, l::Int64)
    return p * (1 - p) / (((1 - 4 * p / 3) ^ 2) * l)
end


## Kimura80 Distance computation internals

@inline function expected_distance(::Type{Kimura80}, a1::AbstractFloat, a2::AbstractFloat)
    return -0.5 * log(a1 * sqrt(a2))
end

@inline function expected_distance(::Type{Kimura80}, a1::AbstractFloat, a2::AbstractFloat,
    gamma::AbstractFloat)
    b = -1 / alpha
    return alpha * ((a1 ^ b) + 0.5 * (a2 ^ b) - 1.5) / 2
end

macro Kimura80var()
    :(return ((c1 * c1 * P + c3 * c3 * Q) - ((c1 * P + c3 * Q) ^ 2)) / L)
end

@inline function variance(::Type{Kimura80}, P::AbstractFloat, Q::AbstractFloat, L::Int, a1::AbstractFloat,
    a2::AbstractFloat)
    c1 = 1 / a1
    c2 = 1 / a2
    c3 = (c1 + c2) / 2
    @Kimura80var
end

@inline function variance(::Type{Kimura80}, P::AbstractFloat, Q::AbstractFloat, L::Int, a1::AbstractFloat,
    a2::AbstractFloat, gamma::AbstractFloat)
    b = -(1 / gamma + 1)
    c1 = a1 ^ b
    c2 = a2 ^ b
    c3 = (c1 + c2) / 2
    @Kimura80var
end


# Distance computation methods
# ----------------------------

"""
Compute the pairwise genetic distances for a set of aligned nucleotide sequences.

The distance measure to compute is determined by the type provided as the first
parameter. The second parameter provides the set of nucleotide sequences.
"""
function distance end

"""
    distance{T<:MutationType,A<:NucleotideAlphabet}(::Type{Count{T}}, seqs::Vector{BioSequence{A}})

Compute the number of mutations of type `T` between a set of sequences in a
pairwise manner.

This method of distance returns a tuple of the number of mutations of type `T`
between sequences and the number of valid (i.e. non-ambiguous sites) counted by
the function.
"""
function distance{T<:MutationType,A<:NucleotideAlphabet}(::Type{Count{T}}, seqs::Vector{BioSequence{A}})
    return count_mutations(T, seqs)
end

"""
    distance{T<:MutationType,N<:Nucleotide}(::Type{Count{T}}, seqs::Matrix{N})

Compute the number of mutations of type `T` between a set of sequences in a
pairwise manner.

This method of distance returns a tuple of the number of mutations of type `T`
between sequences and the number of valid (i.e. non-ambiguous sites) counted by
the function.

**Note: This method assumes that the sequences are stored in the `Matrix{N}`
provided as `seqs` in sequence major order i.e. each column of the matrix is one
complete nucleotide sequence.**
"""
function distance{T<:MutationType,N<:Nucleotide}(::Type{Count{T}}, seqs::Matrix{N})
    return count_mutations(T, seqs)
end

function distance{T<:TsTv,A<:NucleotideAlphabet}(::Type{Count{T}}, seqs::Vector{BioSequence{A}})
    return count_mutations(TransitionMutation, TransversionMutation, seqs)
end

function distance{T<:TsTv,N<:Nucleotide}(::Type{Count{T}}, seqs::Matrix{N})
    return count_mutations(TransitionMutation, TransversionMutation, seqs)
end

"""
    distance{T<:MutationType,A<:NucleotideAlphabet}(::Type{Proportion{T}}, seqs::Vector{BioSequence{A}})

This method of distance returns a tuple of a vector of the p-distances, and a
vector of the number of valid (i.e. non-ambiguous sites) counted by the function.
"""
function distance{T<:MutationType,A<:NucleotideAlphabet}(::Type{Proportion{T}}, seqs::Vector{BioSequence{A}})
    d, l = distance(Count{T}, seqs)
    D = Vector{Float64}(length(d))
    @inbounds for i in 1:length(D)
        D[i] = d[i] / l[i]
    end
    return D, l
end

"""
    distance{T<:MutationType,N<:Nucleotide}(::Type{Proportion{T}}, seqs::Matrix{N})

This method of distance returns a tuple of a vector of the p-distances, and a
vector of the number of valid (i.e. non-ambiguous sites) counted by the function.

**Note: This method assumes that the sequences are stored in the `Matrix{N}`
provided as `seqs` in sequence major order i.e. each column of the matrix is one
complete nucleotide sequence.**
"""
function distance{T<:MutationType,N<:Nucleotide}(::Type{Proportion{T}}, seqs::Matrix{N})
    d, l = distance(Count{T}, seqs)
    D = Vector{Float64}(length(d))
    @inbounds for i in 1:length(D)
        D[i] = d[i] / l[i]
    end
    return D, l
end

"""
    distance{A<:NucleotideAlphabet}(::Type{JukesCantor69}, seqs::Vector{BioSequence{A}})

This method of distance returns a tuple of the expected JukesCantor69 distance
estimate, and the computed variance.
"""
function distance{A<:NucleotideAlphabet}(::Type{JukesCantor69}, seqs::Vector{BioSequence{A}})
    p, l = distance(Proportion{AnyMutation}, seqs)
    D = Vector{Float64}(length(p))
    V = Vector{Float64}(length(p))
    @inbounds for i in 1:length(p)
        D[i] = expected_distance(JukesCantor69, p[i])
        V[i] = variance(JukesCantor69, p[i], l[i])
    end
    return D, V
end

"""
    distance{N<:Nucleotide}(::Type{JukesCantor69}, seqs::Matrix{N})

This method of distance returns a tuple of the expected JukesCantor69 distance
estimate, and the computed variance.

**Note: This method assumes that the sequences are stored in the `Matrix{N}`
provided as `seqs` in sequence major order i.e. each column of the matrix is one
complete nucleotide sequence.**
"""
function distance{N<:Nucleotide}(::Type{JukesCantor69}, seqs::Matrix{N})
    p, l = distance(Proportion{AnyMutation}, seqs)
    D = Vector{Float64}(length(p))
    V = Vector{Float64}(length(p))
    @inbounds for i in 1:length(p)
        D[i] = expected_distance(JukesCantor69, p[i])
        V[i] = variance(JukesCantor69, p[i], l[i])
    end
    return D, V
end


"""
    distance{A<:NucleotideAlphabet}(::Type{Kimura80}, seqs::Vector{BioSequence{A}})

This method of distance returns a tuple of the expected Kimura80 distance
estimate, and the computed variance.
"""
function distance{A<:NucleotideAlphabet}(::Type{Kimura80}, seqs::Vector{BioSequence{A}})
    ns, nv, l = distance(Count{Kimura80}, seqs)
    D = Vector{Float64}(length(ns))
    V = Vector{Float64}(length(ns))
    @inbounds for i in 1:length(ns)
        L = l[i]
        P = ns[i] / L
        Q = nv[i] / L
        a1 = 1 - 2 * P - Q
        a2 = 1 - 2 * Q
        D[i] = expected_distance(Kimura80, a1, a2)
        V[i] = variance(Kimura80, P, Q, L, a1, a2)
    end
    return D, V
end

"""
    distance{N<:Nucleotide}(::Type{Kimura80}, seqs::Matrix{N})

This method of distance returns a tuple of the expected Kimura80 distance
estimate, and the computed variance.

**Note: This method assumes that the sequences are stored in the `Matrix{N}`
provided as `seqs` in sequence major order i.e. each column of the matrix is one
complete nucleotide sequence.**
"""
function distance{N<:Nucleotide}(::Type{Kimura80}, seqs::Matrix{N})
    ns, nv, l = distance(Count{Kimura80}, seqs)
    D = Vector{Float64}(length(ns))
    V = Vector{Float64}(length(ns))
    @inbounds for i in 1:length(ns)
        L = l[i]
        P = ns[i] / L
        Q = nv[i] / L
        a1 = 1 - 2 * P - Q
        a2 = 1 - 2 * Q
        D[i] = expected_distance(Kimura80, a1, a2)
        V[i] = variance(Kimura80, P, Q, L, a1, a2)
    end
    return D, V
end
