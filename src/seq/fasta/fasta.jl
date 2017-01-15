# FASTA
# =====
#
# Reader and writer of the FASTA file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"Metadata for FASTA sequence records containing just a `description` field"
type FASTAMetadata
    description::StringField
end

function FASTAMetadata()
    return FASTAMetadata(StringField())
end

function Base.:(==)(a::FASTAMetadata, b::FASTAMetadata)
    return a.description == b.description
end

function Base.copy(metadata::FASTAMetadata)
    return FASTAMetadata(copy(metadata.description))
end

"""
    FASTASeqRecord(name, seq[, description=""])

Create a sequence record for the FASTA file format.
"""
typealias FASTASeqRecord{S} SeqRecord{S,FASTAMetadata}

function (::Type{FASTASeqRecord})(name::AbstractString,
                                  seq::Sequence,
                                  description::AbstractString="")
    return SeqRecord(name, seq, FASTAMetadata(description))
end

include("fai.jl")
include("reader.jl")
include("writer.jl")
