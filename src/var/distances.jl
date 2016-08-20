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

# Evolutionary distances
abstract EvolutionaryDistances
abstract UncorrectedDist <: EvolutionaryDistances
abstract CorrectedDist <: EvolutionaryDistances
abstract TsTv <: CorrectedDist

immutable Count{T} <: UncorrectedDist end
immutable Proportion{T} <: UncorrectedDist end
immutable JukesCantor69 <: CorrectedDist end
immutable Kimura80 <: TsTv end

typealias JC69 JukesCantor69
typealias K80 Kimura80


# Distance computation internals
# ------------------------------


## K80 Distance computation internals

@inline function expected_distance(::Type{K80}, a1::Float64, a2::Float64)
    return -0.5 * log(a1 * sqrt(a2))
end

@inline function expected_distance(::Type{K80}, a1::Float64, a2::Float64,
    gamma::Float64)
    b = -1 / alpha
    return alpha * ((a1 ^ b) + 0.5 * (a2 ^ b) - 1.5) / 2
end

macro k80var()
    :(return ((c1 * c1 * P + c3 * c3 * Q) - ((c1 * P + c3 * Q) ^ 2)) / L)
end

@inline function variance(::Type{K80}, P::AbstractFloat, Q::AbstractFloat, L::Int, a1::AbstractFloat,
    a2::AbstractFloat)
    c1 = 1 / a1
    c2 = 1 / a2
    c3 = (c1 + c2) / 2
    @k80var
end

@inline function variance(::Type{K80}, P::AbstractFloat, Q::AbstractFloat, L::Int, a1::AbstractFloat,
    a2::AbstractFloat, gamma::AbstractFloat)
    b = -(1 / gamma + 1)
    c1 = a1 ^ b
    c2 = a2 ^ b
    c3 = (c1 + c2) / 2
    @k80var
end


# Distance computation methods
# ----------------------------

"""
    distance{T<:MutationType}(::Type{Count{T}}, a::BioSequence, b::BioSequence)

Compute the genetic distance between two DNA or RNA sequences.

Providing the distance type Count{T} results in a distance which is the
count of the mutations of type T that exist between the two sequences.

For description of all the mutation types, see the Var module documentation.

You may also provide the type TsTv as T, in which case the count returned is a sum
of the number of transitions, and the number of transversions.

This function returns a tuple of the number of mutations, and the number of valid
(i.e. non-ambiguous sites tested by the function).
"""
@inline function distance{T<:MutationType}(::Type{Count{T}}, a::BioSequence, b::BioSequence)
    return count_mutations(a, b, T)
end

@inline function distance{T<:TsTv}(::Type{Count{T}}, a::BioSequence, b::BioSequence)
    return count_mutations(a, b, TransitionMutation, TransversionMutation)
end

"""
    distance{T<:MutationType}(::Type{Proportion{T}}, a::BioSequence, b::BioSequence)

Compute the genetic distance between two DNA or RNA sequences.

Providing the distance type Proportion{T} results in a distance which is the
count of the mutations of type T that exist between the two sequences, divided
by the number of sites examined.

In other words this so called p-distance is simply the proportion of sites
between each pair of sequences, that are mutated (again where T determines
what kind of mutation).

For description of all the mutation types, see the Var module documentation.

You may also provide the type TsTv as T, in which case the count returned is a sum
of the number of transitions, and the number of transversions.

This function returns a tuple of the p-distance, and the number of valid
(i.e. non-ambiguous sites tested by the function).
"""
@inline function distance{T<:MutationType}(::Type{Proportion{T}}, a::BioSequence, b::BioSequence)
    d, l = distance(N_Mutations{T}, a, b)
    return d / l, l
end

"""
    distance(::Type{JukesCantor69}, a::BioSequence, b::BioSequence)

Compute the genetic distance between two DNA or RNA sequences.

Providing the distance type JukesCantor69 (alias JC69) results in a
p-distance (the proportion of non-identical sites between the two sequences),
corrected by the substitution model developed by Jukes and Cantor in 1969.

It assumes that all substitutions (i.e. a change of a base by another one) have
the same probability.
This probability is the same for all sites along the DNA sequence.

This function returns a tuple of the expected K80 distance estimate, and the
computed variance.
"""
function distance(::Type{JC69}, a::BioSequence, b::BioSequence)
    p, l = distance(Proportion{DifferentMutation}, a, b)
    @assert 0.0 <= p <= 0.75 throw(DomainError())
    D = -0.75 * log(1 - 4 * p / 3)
    V = p * (1 - p) / (((1 - 4 * p / 3) ^ 2) * l)
    return D, V
end


"""
    distance(::Type{Kimura80}, a::BioSequence, b::BioSequence)

Compute the genetic distance between two DNA or RNA sequences.

Providing the distance type Kimura80 (alias K80) computes a
distance using the substitution model developed by Kimura in 1980.
It is somtimes called Kimura's 2 parameter distance.

The model makes the same assumptions as Jukes and Cantor's model, but
two-kinds of mutation are considered: Transitions and Transversions.
These two mutations can occur with different probabilities in this model.
Both transition and transversion rates are the same for all sites along the DNA
sequence.

This function returns a tuple of the expected K80 distance estimate, and the
computed variance.
"""
function distance(::Type{K80}, a::BioSequence, b::BioSequence)
    ns, nv, l = distance(Count{K80}, a, b)
    P = ns / l
    Q = nv / l
    a1 = 1 - 2 * P - Q
    a2 = 1 - 2 * Q
    D = expected_distance(K80, a1, a2)
    V = variance(K80, P, Q, l, a1, a2)
    return D, V
end
