# The Bio.Var module
# ==================
#
# Types and methods for analysing biological variation.
#
# Part of the Bio.Var module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Var
using Bio.Seq

export

    # Mutation counting
    MutationType,
    DifferentMutation,
    TransitionMutation,
    TransversionMutation,

    count_mutations,

    # Genetic and Evolutionary distances
    EvolutionaryDistance,
    Count,
    Proportion,
    JukesCantor69,
    Kimura80,

    distance




include("mutation_counting.jl")
include("distances.jl")

end # module Var
