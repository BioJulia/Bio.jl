module Intervals

export Strand,
    Interval,
    IntervalCollection,
    IntervalStream,
    STRAND_NA,
    STRAND_POS,
    STRAND_NEG,
    STRAND_BOTH,
    coverage,
    isoverlapping,
    BED,
    BEDMetadata,
    BEDInterval,
    BigBed,
    BigWig

using Base.Intrinsics,
    Bio.Ragel,
    Bio.StringFields,
    BufferedStreams,
    IntervalTrees,
    Switch,
    Libz,
    Colors

using Base.Collections:
    heappush!,
    heappop!

using Bio:
    AbstractParser,
    FileFormat

import Iterators

import Base:
    show,
    isless,
    push!,
    shift!,
    intersect,
    start,
    next,
    done,
    length,
    convert,
    read,
    read!,
    write,
    getindex,
    get,
    isempty,
    endof,
    ==,
    reverse!,
    open,
    eltype,
    copy

import IntervalTrees:
    first,
    last


include("interval.jl")
include("stream_buffer.jl")
include("intervalcollection.jl")
include("intervalstream.jl")

# Parsing file types
include("bed.jl")
include("bigbed.jl")

end # module Intervals
