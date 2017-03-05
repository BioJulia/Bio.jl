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
import Automa
import Automa.RegExp: @re_str
import BGZFStreams: BGZFStream
# TODO: Needs this branch: https://github.com/BioJulia/BufferedStreams.jl/pull/33
import BufferedStreams: BufferedStreams, BufferedInputStream
importall Bio

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

    distance,

    VCFMetaInfo,
    VCFHeader,
    VCFRecord,
    VCFReader,
    VCFWriter,
    BCFRecord,
    BCFReader,
    BCFWriter,
    header,
    metainfotag,
    metainfoval,
    isfilled,
    chromosome,
    leftposition,
    quality,
    filter_,
    identifier,
    reference,
    alternate,
    information,
    infokeys,
    format,
    genotype

# Bio.@reexport import Bio: isfilled, leftposition

include("mutation_counting.jl")
include("distances.jl")
include("vcf/vcf.jl")
include("bcf/bcf.jl")

end # module Var
