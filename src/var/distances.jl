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


## Kimura80 Distance computation internals

@inline function expected_distance(::Type{Kimura80}, a1::Float64, a2::Float64)
    return -0.5 * log(a1 * sqrt(a2))
end

@inline function expected_distance(::Type{Kimura80}, a1::Float64, a2::Float64,
    gamma::Float64)
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
    distance(::Type{EvolutionaryDistance}, a::BioSequence, b::BioSequence)

Compute the genetic distance between two nucleotide sequences of equal length.

The distance measure to compute is determined by the type provided as the first
parameter. The second and third parameter provide the two nucleotide sequences.
"""
function distance end

"""
    distance{T<:MutationType}(::Type{Count{T}}, a::BioSequence, b::BioSequence)

This method of distance returns a tuple of the number of mutations of type T
between two sequences and the number of valid
(i.e. non-ambiguous sites) counted by the function.
"""
@inline function distance{T<:MutationType}(::Type{Count{T}}, a::BioSequence, b::BioSequence)
    return count_mutations(a, b, T)
end

@inline function distance{T<:TsTv}(::Type{Count{T}}, a::BioSequence, b::BioSequence)
    return count_mutations(a, b, TransitionMutation, TransversionMutation)
end

"""
    distance{T<:MutationType}(::Type{Proportion{T}}, a::BioSequence, b::BioSequence)

This method of distance returns a tuple of the p-distance, and the number of valid
(i.e. non-ambiguous sites) counted by the function.
"""
@inline function distance{T<:MutationType}(::Type{Proportion{T}}, a::BioSequence, b::BioSequence)
    d, l = distance(Count{T}, a, b)
    return d / l, l
end

"""
    distance(::Type{JukesCantor69}, a::BioSequence, b::BioSequence)
    distance(::Type{JukesCantor69}, a::BioSequence, b::BioSequence)

This method of distance returns a tuple of the expected JukesCantor69 distance
estimate, and the computed variance.
"""
function distance(::Type{JukesCantor69}, a::BioSequence, b::BioSequence)
    p, l = distance(Proportion{DifferentMutation}, a, b)
    @assert 0.0 <= p <= 0.75 throw(DomainError())
    D = -0.75 * log(1 - 4 * p / 3)
    V = p * (1 - p) / (((1 - 4 * p / 3) ^ 2) * l)
    return D, V
end


"""
    distance(::Type{Kimura80}, a::BioSequence, b::BioSequence)
    distance(::Type{Kimura80}, a::BioSequence, b::BioSequence)

This method of distance returns a tuple of the expected Kimura80 distance
estimate, and the computed variance.
"""
function distance(::Type{Kimura80}, a::BioSequence, b::BioSequence)
    ns, nv, l = distance(Count{Kimura80}, a, b)
    P = ns / l
    Q = nv / l
    a1 = 1 - 2 * P - Q
    a2 = 1 - 2 * Q
    D = expected_distance(Kimura80, a1, a2)
    V = variance(Kimura80, P, Q, l, a1, a2)
    return D, V
end
