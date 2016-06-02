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


# Parser
# ------

"A type encapsulating the current state of a FASTA parser"
type FASTAParser{S<:Sequence} <: AbstractParser
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}

    function FASTAParser(input::BufferedInputStream)
        return new(Ragel.State(fastaparser_start, input),
                   BufferedOutputStream())
    end
end

Base.eltype{S}(::Type{FASTAParser{S}}) = SeqRecord{S,FASTAMetadata}
Base.eof(parser::FASTAParser) = eof(parser.state.stream)

include("fasta-parser.jl")

function Base.open{S}(input::BufferedInputStream, ::Type{FASTA},
                      ::Type{S}=BioSequence)
    return FASTAParser{S}(input)
end

function Base.show{S}(io::IO, seqrec::SeqRecord{S,FASTAMetadata})
    print(io, ">", seqrec.name)
    if !isempty(seqrec.metadata.description)
        print(io, " ", seqrec.metadata.description)
    end
    print(io, '\n')
    print(io, seqrec.seq)
end

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
