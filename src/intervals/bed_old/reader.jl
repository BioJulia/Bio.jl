# Reader
# ======

"""
    BEDReader(input::IO)

Create a data reader of the BED file format.

# Arguments
* `input`: data source
"""
type BEDReader <: Bio.IO.AbstractReader
    state::Ragel.State

    # intermediate values used during parsing
    red::Float32
    green::Float32
    blue::Float32
    block_size_idx::Int
    block_first_idx::Int

    function BEDReader(input::BufferedInputStream)
        return new(Ragel.State(bedparser_start, input), 0, 0, 0, 1, 1)
    end
end

function Bio.IO.stream(reader::BEDReader)
    return reader.state.stream
end

function Intervals.metadatatype(::BEDReader)
    return BEDMetadata
end

function Base.eltype(::Type{BEDReader})
    return BEDInterval
end

function BEDReader(input::IO)
    return BEDReader(BufferedInputStream(input))
end

function IntervalCollection(interval_stream::BEDReader)
    intervals = collect(BEDInterval, interval_stream)
    return IntervalCollection{BEDMetadata}(intervals, true)
end
