
type GFF3Reader <: Bio.IO.AbstractReader
    state::Ragel.State
    version::VersionNumber
    sequence_regions::Vector{Interval{Void}}

    function GFF3Reader(input::BufferedInputStream)
        return new(Ragel.State(gff3parser_start, input), VersionNumber(0), [])
    end
end

function Bio.IO.stream(reader::GFF3Reader)
    return reader.state.stream
end

function Intervals.metadatatype(::GFF3Reader)
    return GFF3Metadata
end

function Base.eltype(::Type{GFF3Reader})
    return GFF3Interval
end

function GFF3Reader(input::IO)
    return GFF3Reader(BufferedInputStream(input))
end

function IntervalCollection(interval_stream::GFF3Reader)
    intervals = collect(GFF3Interval, interval_stream)
    return IntervalCollection{GFF3Metadata}(intervals, true)
end
