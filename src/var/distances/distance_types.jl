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
immutable Proportion{T<:Mutation} <: UncorrectedDistance end

@inline count_site{T<:Mutation}(::Type{Proportion{T}}) = T

@inline function process_count{T<:Mutation}(::Type{Proportion{T}}, counts)
    l = lengthlist(counts)
    D = Vector{Float64}(l)
    list = getlist(counts)
    @inbounds for i in 1:l
        le = list[i]
        D[i] = le[1] / le[2]
    end
    return PairwiseListMatrix(D, false)
end


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
