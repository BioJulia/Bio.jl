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


# Reader
# ------

include("fai.jl")

"A type encapsulating the current state of a FASTA reader"
type FASTAReader{S<:Sequence} <: Bio.IO.AbstractReader
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    index::Nullable{FASTAIndex}

    function FASTAReader(input::BufferedInputStream, index=Nullable())
        return new(Ragel.State(fastaparser_start, input),
                   BufferedOutputStream(), index)
    end
end

function Bio.IO.stream(reader::FASTAReader)
    return reader.state.stream
end

Base.eltype{S}(::Type{FASTAReader{S}}) = FASTASeqRecord{S}

function Base.getindex(reader::FASTAReader, name::AbstractString)
    if isnull(reader.index)
        error("no index")
    end
    seekrecord(reader.state.stream, get(reader.index), name)
    reader.state.cs = fastaparser_start
    return read(reader)
end

include("fasta-parser.jl")


# Writer
# ------

# Serializer for the FASTA file format.
type FASTAWriter{T<:IO} <: Bio.IO.AbstractWriter
    output::T
    # maximum sequence width (no limit when width ≤ 0)
    width::Int
end

function Bio.IO.stream(writer::FASTAWriter)
    return writer.output
end

function Base.write(writer::FASTAWriter, seqrec::FASTASeqRecord)
    output = writer.output
    n = 0

    # header
    n += write(output, '>', seqrec.name)
    if !isempty(seqrec.metadata.description)
        n += write(output, ' ', seqrec.metadata.description)
    end
    n += write(output, '\n')

    # sequence
    w = writer.width
    for x in seqrec.seq
        if writer.width > 0 && w == 0
            n += write(output, '\n')
            w = writer.width
        end
        n += write(output, Char(x))
        w -= 1
    end
    n += write(output, '\n')

    return n
end
