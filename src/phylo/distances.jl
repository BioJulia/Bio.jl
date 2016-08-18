# phylo/distances.jl
# ==================
#
# Types and methods for computing evolutionary and genetic distances.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Types
# -----

# Evolutionary distances
abstract EvoDist
abstract UncorrectedDist <: EvoDist
abstract CorrectedDist <: EvoDist
abstract TsTv <: CorrectedDist

immutable N_Mutations{T} <: UncorrectedDist end
immutable P_Distance{T} <: UncorrectedDist end
immutable JukesCantor69 <: CorrectedDist end
immutable Kimura80 <: TsTv end

typealias Raw N_Mutations{DifferentMutation}
typealias P P_Distance{DifferentMutation}
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

# Distance method that essentially just calls the count_mutations function.
@inline function distance{T<:MutationType}(::Type{N_Mutations{T}}, a::BioSequence, b::BioSequence)
    return count_mutations(a, b, T)
end

@inline function distance{T<:TsTv}(::Type{N_Mutations{T}}, a::BioSequence, b::BioSequence)
    return count_mutations(a, b, TransitionMutation, TransversionMutation)
end

# Method for computing the P distance of any kind of mutation.
@inline function distance{T<:MutationType}(::Type{P_Distance{T}}, a::BioSequence, b::BioSequence)
    d, l = distance(N_Mutations{T}, a, b)
    return d / l
end

function distance(::Type{JC69}, a::BioSequence, b::BioSequence)
    n, l = distance(N_Mutations{DifferentMutation}, a, b)
    p = n / l
    @assert 0.0 <= p <= 0.75 throw(DomainError())
    D = -0.75 * log(1 - 4 * p / 3)
    V = p * (1 - p) / (((1 - 4 * p / 3) ^ 2) * l)
    return D, V
end

function distance(::Type{JC69}, a::BioSequence, b::BioSequence, alpha::Float64)
    n, l = distance(N_Mutations{DifferentMutation}, a, b)
    p = n / l
    @assert 0.0 <= p <= 0.75 throw(DomainError())
    D = 0.75 * alpha * ( (1 - 4 * p / 3) ^ (-1 / alpha) - 1)
    V = p * (1 - p)/(((1 - 4 * p / 3) ^ (-2 / (alpha + 1))) * l)
    return D, V
end

function distance(::Type{K80}, a::BioSequence, b::BioSequence)
    ns, nv, l = distance(N_Mutations{K80}, a, b)
    P = ns / l
    Q = nv / l
    a1 = 1 - 2 * P - Q
    a2 = 1 - 2 * Q
    D = expected_distance(K80, a1, a2)
    V = variance(K80, P, Q, l, a1, a2)
    return D, V
end

function distance(::Type{K80}, a::BioSequence, b::BioSequence, gamma::Float64)
    ns, nv, l = distance(N_Mutations{K80}, a, b)
    P = ns / l
    Q = nv / l
    a1 = 1 - 2 * P - Q
    a2 = 1 - 2 * Q
    D = expected_distance(K80, a1, a2, gamma)
    V = variance(K80, P, Q, l, a1, a2, gamma)
    return D, V
end
