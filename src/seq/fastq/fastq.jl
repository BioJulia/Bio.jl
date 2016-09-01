# FASTQ
# =====
#
# Reader and writer of the FASTQ file format.
#
# Cock, Peter JA, et al. "The Sanger FASTQ file format for sequences with
# quality scores, and the Solexa/Illumina FASTQ variants." Nucleic acids
# research 38.6 (2010): 1767-1771.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
Metadata for FASTQ sequence records containing a `description` field,
and a `quality` string corresponding to the sequence.

Quality scores are stored as integer Phred scores.
"""
type FASTQMetadata
    description::StringField
    quality::Vector{Int8}
end

function FASTQMetadata()
    return FASTQMetadata(StringField(), Int8[])
end

function Base.:(==)(a::FASTQMetadata, b::FASTQMetadata)
    return a.description == b.description && a.quality == b.quality
end

function Base.copy(metadata::FASTQMetadata)
    return FASTQMetadata(copy(metadata.description), copy(metadata.quality))
end

"""
    FASTQSeqRecord(name, seq, quality[, description=""])

Create a sequence record for the FASTQ file format.
"""
typealias FASTQSeqRecord{S} SeqRecord{S,FASTQMetadata}

function (::Type{FASTQSeqRecord})(name::AbstractString,
                                  seq::Sequence,
                                  quality::Vector{Int8},
                                  description::AbstractString="")
    if length(seq) != length(quality)
        throw(ArgumentError("sequence and base quality must be the same length"))
    end
    return SeqRecord(name, seq, FASTQMetadata(description, quality))
end

include("quality.jl")
include("reader.jl")
include("parser.jl")
include("writer.jl")
