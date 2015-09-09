module Intervals

import Base: show, isless, push!, shift!, intersect, start, next, done, length,
             convert, read, read!, write, getindex, get, isempty, endof, ==,
             reverse!

import Bio: FileFormat
import IntervalTrees: first, last
import Iterators, Zlib

#using DataStructures
using Base.Intrinsics, Color, IntervalTrees


export Strand, Interval, IntervalCollection, IntervalStream,
       STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH,
       coverage, isoverlapping, BED, BEDMetadata, BEDInterval, BigBed, BigWig

include("interval.jl")
include("stream_buffer.jl")
include("intervalcollection.jl")
include("intervalstream.jl")

# Parsing file types
include("bed.jl")

# TODO: reorg
include("bigbed.jl")

end
