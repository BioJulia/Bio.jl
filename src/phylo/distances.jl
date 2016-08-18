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

## JC69 Distance computation internals.

@inline function expected_distance(::Type{JC69}, x::Float64)
    return -0.75 * log(1 - 4 * x / 3)
end

@inline function expected_distance(::Type{JC69}, x::Float64,
    gamma::Float64)
    return 0.75 * alpha * ( (1 - 4 * p / 3) ^ (-1 / alpha) - 1)
end

@inline function variance(::Type{JC69}, x::Float64, L::Int)
    return x * (1 - x) / (((1 - 4 * p / 3) ^ 2) * L)
end

@inline function variance(::Type{JC69}, x::Float64, L::Int,
    gamma::Float64)
    return x * (1 - x)/(((1 - 4 * x / 3) ^ (-2 / (alpha + 1))) * L)
end

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
    :(return ((c1 * c1 * P + c3 * c3 * Q) - ((c1 * P + c3 * q) ^ 2)) / L)
end

@inline function variance(::Type{K80}, P::Int, Q::Int, L::Int, a1::Float64,
    a2::Float64)
    c1 = 1 / a1
    c2 = 1 / a2
    c3 = (c1 + c2) / 2
    @k80var
end

@inline function variance(::Type{K80}, P::Int, Q::Int, L::Int, a1::Float64,
    a2::Float64, gamma::Float64)
    b = -(1 / alpha + 1)
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

@inline function distance(::Type{N_Mutations{TsTv}}, a::BioSequence, b::BioSequence)
    return count_mutations(a, b, TransitionMutation, TransversionMutation)
end

# Method for computing the P distance of any kind of mutation.
@inline function distance{T<:MutationType}(::Type{P_Distance{T}}, a::BioSequence, b::BioSequence)
    d, l = distance(N_Mutations{T}, a, b)
    return d / l
end

# Method to compute distance corrected by JukesCantor69 substitution model.
function distance(::Type{JC69}, a::BioSequence, b::BioSequence)
    p = distance(P_Distance{DifferentMutation}, a, b)
    D = expected_distance(JukesCantor69, p)
    V = variance(JukesCantor69, p, l)
    return D, V
end

function distance(::Type{JC69}, a::BioSequence, b::BioSequence, gamma::Float64)
    p = distance(P_Distance{DifferentMutation}, a, b)
    D = expected_distance(JukesCantor69, p, gamma)
    V = variance(JukesCantor69, p, l, gamma)
    return D, V
end














function distance(::Type{Kimura80}, a::BioSequence, b::BioSequence)
    nd, ns, nv, l = countTsTv(a, b, model)
    P = ns / L
    Q = nv / L
    a1 = 1 - 2 * P - Q
    a2 = 1 - 2 * Q
    D = model_correction(model, a1, a2)
    V = variance(model, P, Q, L, a1, a2)
    return D, V
end

function distance(a::BioSequence, b::BioSequence, t::Type{Kimura80}, gamma::Float64)
    nd, ns, nv, l = countTsTv(a, b, model)
    P = ns / L
    Q = nv / L
    a1 = 1 - 2 * P - Q
    a2 = 1 - 2 * Q
    D = model_correction(model, a1, a2, gamma)
    V = variance(model, P, Q, L, a1, a2, gamma)
    return D, V
end
