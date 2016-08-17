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

immutable N_Mutations{T<:MutationType} <: UncorrectedDist end
immutable P_Distance{T<:MutationType} <: UncorrectedDist end
immutable JukesCantor69 <: CorrectedDist end
immutable Kimura80 <: TsTv end


# Distance computation methods
# ----------------------------

# Indicate inline on these functions as they are called by other functions repeatedly in a loop.

## Number of mutations.

# Conveinience method that essentially just calls the count_mutations function.
@inline function distance{T<:MutationType}(a::BioSequence, b::BioSequence, t::Type{N_Mutations{T}})
    d, l = count_mutations(a, b, T)
    return d, l
end

# Method for computing the P distance of any kind of mutation.
@inline function distance{T<:MutationType}(a::BioSequence, b::BioSequence, t::Type{P_Distance{T}})
    d, l = count_mutations(a, b, count_type(t))
    return d / l
end

# Method to compute distance corrected by JukesCantor69 substitution model.
@inline function distance(a::BioSequence, b::BioSequence, t::Type{JukesCantor69})
    d, l = count_mutations(a, b, DifferentMutation)
    p = d / l
    D = expected_distance(p, t)
    V = variance(p, l, t)
    return D, V
end







## JC69 Distance computation

function expected_distance(x::Float64, t::Type{JukesCantor69})
    return -0.75 * log(1 - 4 * x / 3)
end

function expected_distance(x::Float64, gamma::Float64, t::Type{JukesCantor69})
    return 0.75 * alpha * ( (1 - 4 * p / 3) ^ (-1 / alpha) - 1)
end

function variance(x::Float64, L::Int, t::Type{JukesCantor69})
    return x * (1 - x) / (((1 - 4 * p / 3) ^ 2) * L)
end

function variance(x::Float64, L::Int, gamma::Float64, t::Type{JukesCantor69})
    return x * (1 - x)/(((1 - 4 * x / 3) ^ (-2 / (alpha + 1))) * L)
end



function distance(a::BioSequence, b::BioSequence, gamma::Float64, t::Type{JukesCantor69})
    d, l = count_differences(a, b, model)
    p = d / l
    D = expected_distance(p, gamma, model)
    V = variance(p, l, gamma, model)
    return D, V
end

function distance_pairdel(a::BioSequence, b::BioSequence, t::Type{JukesCantor69})
    d, l = count_differences_pairdel(a, b, model)
    p = d / l
    D = expected_distance(p, model)
    V = variance(p, l, model)
    return D, V
end

function distance_pairdel(a::BioSequence, b::BioSequence, t::Type{JukesCantor69}, gamma::Float64)
    d, l = count_differences_pairdel(a, b, model)
    p = d / l
    D = model_correction(p, gamma, model)
    V = variance(p, l, gamma, model)
    return D, V
end



# K80 Distance computation
# ------------------------

function model_correction(t::Type{Kimura80}, a1::Float64, a2::Float64)
    return -0.5 * log(a1 * sqrt(a2))
end

function model_correction(t::Type{Kimura80}, a1::Float64, a2::Float64, gamma::Float64)
    b = -1 / alpha
    return alpha * ((a1 ^ b) + 0.5 * (a2 ^ b) - 1.5) / 2
end

macro k80var()
    :(return ((c1 * c1 * P + c3 * c3 * Q) - ((c1 * P + c3 * q) ^ 2)) / L)
end

function variance(t::Type{Kimura80}, P::Int, Q::Int, L::Int, a1::Float64, a2::Float64)
    c1 = 1 / a1
    c2 = 1 / a2
    c3 = (c1 + c2) / 2
    @k80var
end

function variance(t::Type{Kimura80}, P::Int, Q::Int, L::Int, a1::Float64, a2::Float64, gamma::Float64)
    b = -(1 / alpha + 1)
    c1 = a1 ^ b
    c2 = a2 ^ b
    c3 = (c1 + c2) / 2
    @k80var
end

function distance(a::BioSequence, b::BioSequence, t::Type{Kimura80})
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
