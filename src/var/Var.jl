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

    ## Mutation Counting

    ### Mutation count types.
    DifferentMutation,
    TransitionMutation,
    TransversionMutation,

    ### Checking mutation types of nucleotides.

    count_mutations

include("mutation_counting.jl")

end # module Var
