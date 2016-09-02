# FASTQ Reader
# ============

"A type encapsulating the current state of a FASTQ reader"
type FASTQReader{S<:Sequence} <: Bio.IO.AbstractReader
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    qualbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    name2buf::StringField
    desc2buf::StringField
    qualcount::Int
    quality_encoding::QualityEncoding
    fillN::DNANucleotide

    function FASTQReader(input::BufferedInputStream, quality_encoding::QualityEncoding, fillN::DNANucleotide)
        return new(Ragel.State(fastqparser_start, input),
                   BufferedOutputStream(), BufferedOutputStream(),
                   StringField(), StringField(), 0, quality_encoding, fillN)
    end
end

function Bio.IO.stream(reader::FASTQReader)
    return reader.state.stream
end

function Base.eltype{S}(::Type{FASTQReader{S}})
    return FASTQSeqRecord{S}
end

function FASTQReader(input::IO; quality_encoding::Symbol=:sanger, fillN::DNANucleotide=DNA_N)
    return FASTQReader{DNASequence}(BufferedInputStream(input), sym2qualenc(quality_encoding), fillN)
end

function (::Type{FASTQReader{S}}){S<:Sequence}(input::IO;
                                               quality_encoding::Symbol=:sanger,
                                               fillN::DNANucleotide=DNA_N)
    return FASTQReader{S}(BufferedInputStream(input), sym2qualenc(quality_encoding), fillN)
end

# Map the quality encoding name to a corresponding QualityEncoding object.
function sym2qualenc(name::Symbol)
    if name == :sanger
        return SANGER_QUAL_ENCODING
    elseif name == :solexa
        return SOLEXA_QUAL_ENCODING
    elseif name == :illumina13
        return ILLUMINA13_QUAL_ENCODING
    elseif name == :illumina15
        return ILLUMINA15_QUAL_ENCODING
    elseif name == :illumina18 || name == :illumina
        return ILLUMINA18_QUAL_ENCODING
    end
    error("quality encoding ':$(name)' is not supported")
end
