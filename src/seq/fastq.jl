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


# Reader
# ------

"A type encapsulating the current state of a FASTQ reader"
type FASTQReader{S<:Sequence} <: Bio.IO.AbstractReader
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    qualbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    name2buf::StringField
    desc2buf::StringField
    qualcount::Int
    quality_encodings::QualityEncoding

    function FASTQReader(input::BufferedInputStream,
                         quality_encodings::QualityEncoding)
        if quality_encodings == EMPTY_QUAL_ENCODING
            error("The `quality_encodings` argument is required when parsing FASTQ.")
        elseif count_ones(convert(UInt16, quality_encodings)) > 1
            error("The `quality_encodings` argument must specify exactly one encoding.")
        elseif count_ones(convert(UInt16, quality_encodings & ALL_QUAL_ENCODINGS)) != 1
            error("Unknown quality encoding.")
        end
        return new(Ragel.State(fastqparser_start, input),
                   BufferedOutputStream(), BufferedOutputStream(),
                   StringField(), StringField(), 0, quality_encodings)
    end
end

function Bio.IO.stream(reader::FASTQReader)
    return reader.state.stream
end

Base.eltype{S}(::Type{FASTQReader{S}}) = FASTQSeqRecord{S}

include("fastq-parser.jl")


# Writer
# ------

type FASTQWriter{T<:IO} <: Bio.IO.AbstractWriter
    output::T

    # write sequence name and description header at the third line
    quality_header::Bool

    # encoding offset of base qualities (ASCII code = Phred quality + offset)
    ascii_offset::Int
end

function Bio.IO.stream(writer::FASTQWriter)
    return writer.output
end

function Base.write(writer::FASTQWriter, seqrec::FASTQSeqRecord)
    if writer.ascii_offset == typemin(Int)
        # infer quality encoding based on data
        if !isempty(seqrec.metadata.quality) && minimum(seqrec.metadata.quality) < 0
            writer.ascii_offset = 64  # solexa
        else
            writer.ascii_offset = 33  # others
        end
    end

    output = writer.output
    n = 0

    # header
    n += write(output, '@', seqrec.name)
    if !isempty(seqrec.metadata.description)
        n += write(output, ' ', seqrec.metadata.description)
    end
    n += write(output, '\n')

    # sequence
    for x in seqrec.seq
        n += write(output, Char(x))
    end
    n += write(output, '\n')

    # (optional) identifier and description
    n += write(output, '+')
    if writer.quality_header
        n += write(output, seqrec.name)
        if !isempty(seqrec.metadata.description)
            n += write(output, ' ', seqrec.metadata.description)
        end
    end
    n += write(output, '\n')

    for q in seqrec.metadata.quality
        n += write(output, Char(q + writer.ascii_offset))
    end
    n += write(output, '\n')

    return n
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
