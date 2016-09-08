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
using Bio.Seq: BitIndex, bitindex, index, offset, mask, bitsof, seq_data_len

export

    # Site Counting
    ## Site types
    SiteCase,
    Conserved,
    Mutated,
    Transition,
    Transversion,
    Gap,
    Ambiguous,
    Pairdel



    # Identifying and counting mutations
    #count_mutations,
    #is_mutation,
    #flagmutations,

    # Genetic and Evolutionary distances
    #EvolutionaryDistance,
    #Count,
    #Proportion,
    #JukesCantor69,
    #Kimura80,

    #distance



include("site_counting/site_counting.jl")


end # module Var
