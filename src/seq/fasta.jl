# FASTA
# =====
#
# Reader and writer of the FASTA file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable FASTA <: FileFormat end

"Metadata for FASTA sequence records containing just a `description` field"
type FASTAMetadata
    description::StringField
end

function FASTAMetadata()
    return FASTAMetadata(StringField())
end

@compat function Base.:(==)(a::FASTAMetadata, b::FASTAMetadata)
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

@compat function (::Type{FASTASeqRecord})(name::AbstractString,
                                          seq::Sequence,
                                          description::AbstractString="")
    return SeqRecord(name, seq, FASTAMetadata(description))
end

function Base.open(filepath::AbstractString, ::Type{FASTA})
    input = BufferedInputStream(open(filepath))
    indexpath = filepath * ".fai"
    if isfile(indexpath)
        return FASTAParser{BioSequence}(input, FASTAIndex(indexpath))
    else
        return FASTAParser{BioSequence}(input)
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
    return FASTAParser{S}(input)
end


# Parser
# ------

include("fai.jl")

"A type encapsulating the current state of a FASTA parser"
type FASTAParser{S<:Sequence} <: AbstractParser
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    index::Nullable{FASTAIndex}

    function FASTAParser(input::BufferedInputStream, index=Nullable())
        return new(Ragel.State(fastaparser_start, input),
                   BufferedOutputStream(), index)
    end
end

Base.eltype{S}(::Type{FASTAParser{S}}) = FASTASeqRecord{S}
Base.eof(parser::FASTAParser) = eof(parser.state.stream)

function Base.getindex(parser::FASTAParser, name::AbstractString)
    if isnull(parser.index)
        error("no index")
    end
    seekrecord(parser.state.stream, get(parser.index), name)
    parser.state.cs = fastaparser_start
    return read(parser)
end

include("fasta-parser.jl")


# Writer
# ------

# Serializer for the FASTA file format.
type FASTAWriter{T<:IO} <: AbstractWriter
    output::T
    # maximum sequence width (no limit when width ≤ 0)
    width::Int
end

function Base.flush(writer::FASTAWriter)
    # TODO: This can be removed on Julia v0.5
    # (because flush will be defined for IOBuffer).
    if applicable(flush, writer.output)
        flush(writer.output)
    end
end

Base.close(writer::FASTAWriter) = close(writer.output)

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

function Base.show(io::IO, seqrec::FASTASeqRecord)
    print(io, ">", seqrec.name)
    if !isempty(seqrec.metadata.description)
        print(io, " ", seqrec.metadata.description)
    end
    print(io, '\n')
    print(io, seqrec.seq)
end

# This function is almost deprecated in favor of FASTAWriter.
function Base.write{T}(io::IO, seqrec::SeqRecord{T,FASTAMetadata})
    write(io, ">", seqrec.name)
    if !isempty(seqrec.metadata.description)
        write(io, " ", seqrec.metadata.description)
    end
    write(io, "\n")
    maxchars = 79
    counter = 1
    len = length(seqrec.seq)
    for nt in seqrec.seq
        print(io, nt)
        if counter % maxchars == 0 && counter < len
            write(io, "\n")
        end
        counter += 1
    end
    write(io, "\n")
end
