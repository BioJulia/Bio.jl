# FASTQ
# =====
#
# Reader and writer of the FASTQ file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable FASTQ <: FileFormat end

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

@compat function Base.:(==)(a::FASTQMetadata, b::FASTQMetadata)
    return a.description == b.description && a.quality == b.quality
end

function Base.copy(metadata::FASTQMetadata)
    return FASTQMetadata(copy(metadata.description), copy(metadata.quality))
end


# Parser
# ------

"A type encapsulating the current state of a FASTQ parser"
type FASTQParser{S<:Sequence} <: AbstractParser
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    qualbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    name2buf::StringField
    desc2buf::StringField
    qualcount::Int
    quality_encodings::QualityEncoding

    function FASTQParser(input::BufferedInputStream,
                         quality_encodings::QualityEncoding)
        return new(Ragel.State(fastqparser_start, input),
                   BufferedOutputStream(), BufferedOutputStream(),
                   StringField(), StringField(), 0, quality_encodings)
    end
end

Base.eltype{S}(::Type{FASTQParser{S}}) = SeqRecord{S,FASTQMetadata}
Base.eof(parser::FASTQParser) = eof(parser.state.stream)

include("fastq-parser.jl")

function Base.open{S}(input::BufferedInputStream, ::Type{FASTQ},
                      ::Type{S}=DNASequence;
                      qualenc::QualityEncoding=EMPTY_QUAL_ENCODING)
    return FASTQParser{S}(input, qualenc)
end

"A `SeqRecord` for FASTQ sequences"
typealias FASTQSeqRecord DNASeqRecord{FASTQMetadata}

function FASTQSeqRecord()
    return FASTQSeqRecord(StringField(), DNASequence(), FASTQMetadata())
end

function Base.show(io::IO, seqrec::FASTQSeqRecord)
    write(io, "@", seqrec.name, " ", seqrec.metadata.description, "\n")
    for c in seqrec.seq
        print(io, c)
    end
    write(io, '\n')
    # print quality scores as a unicode bar chart
    for q in seqrec.metadata.quality
        if q <= 0
            write(io, '▁')
        elseif q <= 6
            write(io, '▂')
        elseif q <= 12
            write(io, '▃')
        elseif q <= 18
            write(io, '▄')
        elseif q <= 24
            write(io, '▅')
        elseif q <= 30
            write(io, '▆')
        elseif q <= 36
            write(io, '▇')
        else
            write(io, '█')
        end
    end
    write(io, '\n')
end

function Base.write(io::IO, seqrec::FASTQSeqRecord;
                    offset::Integer=-1, qualheader::Bool=false)

    # choose offset automatically
    if offset < 0
        if !isempty(seqrec.metadata.quality) && minimum(seqrec.metadata.quality) < 0
            offset = 64 # solexa quality offset
        else
            offset = 33  # sanger
        end
    end

    write(io, "@", seqrec.name)
    if !isempty(seqrec.metadata.description)
        write(io, " ", seqrec.metadata.description)
    end
    write(io, "\n")

    for c in seqrec.seq
        print(io, c)
    end
    write(io, "\n")

    write(io, "+")
    if qualheader
        write(io, seqrec.name)
        if !isempty(seqrec.metadata.description)
            write(io, " ", seqrec.metadata.description)
        end
    end
    write(io, "\n")

    for q in seqrec.metadata.quality
        write(io, Char(q + offset))
    end
    write(io, "\n")
end
