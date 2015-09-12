
module Intervals

import Base: show, isless, push!, shift!, intersect, start, next, done, length,
             convert, read, read!, write, getindex, get, isempty, endof, ==,
             reverse!, open, eltype, copy

import Iterators
using Base.Intrinsics
using BufferedStreams
using Colors
using Compat
using IntervalTrees
using Switch
import IntervalTrees: first, last

using Bio: AbstractParser, FileFormat
using Bio.StringFields
using Bio.Ragel

export Strand, Interval, IntervalCollection, IntervalStream,
       STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH,
       coverage, isoverlapping, BED, BEDMetadata, BEDInterval, BigBed, BigWig

include("interval.jl")
include("stream_buffer.jl")
include("intervalcollection.jl")
include("intervalstream.jl")

# Parsing file types
include("bed.jl")
include("bigbed.jl")

end
