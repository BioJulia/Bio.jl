# Bio.Intervals
# =============
#
# Module for genomic intervals.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Intervals

export
    Strand,
    Interval,
    IntervalCollection,
    IntervalStream,
    STRAND_NA,
    STRAND_POS,
    STRAND_NEG,
    STRAND_BOTH,
    coverage,
    isoverlapping,
    seqname,
    metadata,
    strand,
    BED,
    BEDMetadata,
    BEDInterval,
    BigBed,
    BigWig,
    Tabix,
    overlapchunks

import ..Ragel: tryread!
export tryread!

import Bio
using
    BufferedStreams,
    Colors,
    IntervalTrees,
    Libz,
    Bio.Ragel,
    Bio.StringFields,
    Bio.Seq

using Base.Collections:
    heappush!,
    heappop!

import Iterators

include("strand.jl")
include("interval.jl")
include("stream_buffer.jl")
include("intervalcollection.jl")
include("intervalstream.jl")
include("index/index.jl")

# Parsing file types
include("bed/bed.jl")
include("bigbed/bigbed.jl")

end # module Intervals
