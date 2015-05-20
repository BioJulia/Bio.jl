
module Intervals

import Base: show, isless, push!, shift!, intersect, start, next, done, length,
             isempty, endof, getindex
import Iterators
using Base.Intrinsics
#using DataStructures
using Docile
using Docile: @doc, @doc_str
using IntervalTrees
import IntervalTrees: first, last

export Strand, Interval, IntervalCollection,
       STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH,
       isoverlapping

include("interval.jl")
include("stream_buffer.jl")
include("intervalcollection.jl")
include("intervalstream.jl")


end

