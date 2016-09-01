# FASTQ
# =====
#
# Reader and writer of the FASTQ file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable FASTQ <: Bio.IO.FileFormat end

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

function Base.open(
        filepath::AbstractString,
        mode::AbstractString,
        ::Type{FASTQ};
        # reader options
        quality_encoding::QualityEncoding=EMPTY_QUAL_ENCODING,
        # writer options
        quality_header::Bool=false,
        ascii_offset::Integer=typemin(Int))
    io = open(filepath, mode)
    if mode[1] == 'r'
        return open(BufferedInputStream(io), FASTQ, quality_encoding)
    elseif mode[1] ∈ ('w', 'a')
        return FASTQWriter(io, quality_header, ascii_offset)
    end
    error("invalid open mode")
end

function Base.open{S}(input::BufferedInputStream, ::Type{FASTQ},
                      quality_encoding::QualityEncoding,
                      ::Type{S}=DNASequence;
                      # TODO: remove this option after v0.2
                      qualenc=quality_encoding)
    return FASTQReader{S}(input, quality_encoding)
end
