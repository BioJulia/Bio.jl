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

function Base.eltype{S}(::Type{FASTQReader{S}})
    return FASTQSeqRecord{S}
end

function FASTQReader(input::IO, quality_encoding::QualityEncoding)
    return FASTQReader{DNASequence}(BufferedInputStream(input), quality_encoding)
end

function (::Type{FASTQReader{S}}){S<:Sequence}(input::IO, quality_encoding::QualityEncoding)
    return FASTQReader{S}(BufferedInputStream(input), quality_encoding)
end
