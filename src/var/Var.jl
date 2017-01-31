# The Bio.Var module
# ==================
#
# Types and methods for analysing biological variation.
#
# Part of the Bio.Var module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Var

importall Bio
using Bio.Seq

export

    # Mutation types
    MutationType,
    AnyMutation,
    TransitionMutation,
    TransversionMutation,

    # Identifying and counting mutations
    count_mutations,
    is_mutation,
    flagmutations,

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
