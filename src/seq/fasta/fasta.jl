# FASTA
# =====
#
# Reader and writer of the FASTA file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable FASTA <: Bio.IO.FileFormat end

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

function Base.open(filepath::AbstractString, ::Type{FASTA})
    input = BufferedInputStream(open(filepath))
    indexpath = filepath * ".fai"
    if isfile(indexpath)
        return FASTAReader{BioSequence}(input, FASTAIndex(indexpath))
    else
        return FASTAReader{BioSequence}(input)
    end
end

function Base.open(filepath::AbstractString, mode::AbstractString, ::Type{FASTA};
                   width::Integer=60)
    io = open(filepath, mode)
    if mode[1] == 'r'
        return open(BufferedInputStream(io), FASTA)
    elseif mode[1] ∈ ('w', 'a')
        return FASTAWriter(io, width)
    end
    error("invalid open mode")
end

function Base.open{S}(input::BufferedInputStream, ::Type{FASTA},
                      ::Type{S}=BioSequence)
    return FASTAReader{S}(input)
end

include("fai.jl")
include("reader.jl")
include("parser.jl")
include("writer.jl")
