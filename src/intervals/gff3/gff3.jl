
# Generic Feature Format Version 3 (GFF3)
# =======================================
#
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

type GFF3Metadata
    source::StringField
    # use "kind" instead of "type" since "type" is a keyword of Julia
    kind::StringField
    score::Nullable{Float64}
    phase::Nullable{Int}
    attributes::StringField
end

function GFF3Metadata()
    return GFF3Metadata("", "", NaN, 0, "")
end

function Base.copy(metadata::GFF3Metadata)
    return GFF3Metadata(
        copy(metadata.source),
        copy(metadata.kind),
        metadata.score,
        metadata.phase,
        copy(metadata.attributes)
    )
end

"An `Interval` with associated metadata from a GFF3 file"
typealias GFF3Interval Interval{GFF3Metadata}

include("reader.jl")
include("parser.jl")
