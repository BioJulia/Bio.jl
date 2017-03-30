# The Bio.Var module
# ==================
#
# Types and methods for analysing biological variation.
#
# Part of the Bio.Var module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Var

import Bio.Seq:
    Site,
    issite,
    Certain,
    Match,
    Mismatch,
    NucleicAcid,
    MinHashSketch,
    start_count,
    update_count,
    ispurine,
    ispyrimidine
import PairwiseListMatrices: PairwiseListMatrix
import Bio.Exceptions: MissingFieldException, missingerror
import Bio.Windows: eachwindow, EachWindowIterator, SeqWinItr
import Automa
import Automa.RegExp: @re_str
import BGZFStreams: BGZFStream
# TODO: Needs this branch: https://github.com/BioJulia/BufferedStreams.jl/pull/33
import BufferedStreams: BufferedStreams, BufferedInputStream
import IntervalTrees: Interval, IntervalValue
importall Bio

export
    # Site types
    Mutation,
    Conserved,
    Mutated,
    Transition,
    Transversion,

    # VCF and BCF
    VCF,
    BCF,
    header,
    metainfotag,
    metainfoval,
    isfilled,
    chromosome,
    haschromosome,
    leftposition,
    hasleftposition,
    quality,
    hasquality,
    filter_,
    hasfilter_,
    identifier,
    hasidentifier,
    reference,
    hasreference,
    alternate,
    hasalternate,
    information,
    hasinformation,
    infokeys,
    format,
    hasformat,
    genotype,

    MissingFieldException,
    mashdistance

# Bio.@reexport import Bio: isfilled, leftposition

include("site_counting/site_counting.jl")
#include("distances.jl")
include("vcf/vcf.jl")
include("bcf/bcf.jl")
include("mash.jl")

end # module Var
