
type GFF3Reader <: Bio.IO.AbstractReader
    state::Ragel.State
    version::VersionNumber
    sequence_regions::Vector{Interval{Void}}
    key::StringField
    entry_seen::Bool
    fasta_seen::Bool

    function GFF3Reader(input::BufferedInputStream)
        return new(Ragel.State(gff3parser_start, input), VersionNumber(0), [],
                   StringField(), false, false)
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

function hasfasta(reader::GFF3Reader)
    if eof(reader)
        return reader.fasta_seen
    else
        error("GFF3 file must be read until the end before any FASTA sequences can be accessed")
    end
end

function getfasta(reader::GFF3Reader)
    if !hasfasta(reader)
        error("GFF3 file has no FASTA data ")
    end
    return FASTAReader(reader.state.stream)
end
