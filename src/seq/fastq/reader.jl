# FASTQ Reader
# ============

"""
    FASTQReader(input::IO, quality_encoding=:sanger, fill_ambiguous=nothing)

Create a data reader of the FASTQ file format.

# Arguments
* `input`: data source
* `quality_encoding=:sanger`: encoding of base qualities; see the following table for available values
* `fill_ambiguous=nothing`: fill ambiguous nucleotides with the given nucleotide

| Quality encoding | Symbol           | ASCII offset | Quality range |
|:---------------- |:---------------- |:------------:|:-------------:| 
|  Sanger          | `:sanger`        |     +33      |      0-93     |
|  Solexa          | `:solexa`        |     +64      |     -5-62     |
|  Illumina 1.3+   | `:illumina13`    |     +64      |      0-62     |
|  Illumina 1.5+   | `:illumina15`    |     +64      |      2-62     |
|  Illumina 1.8+   | `:illumina18`    |     +33      |      0-93     |
"""
type FASTQReader{S<:Sequence} <: Bio.IO.AbstractReader
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    qualbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    name2buf::StringField
    desc2buf::StringField
    qualcount::Int
    quality_encoding::QualityEncoding
    fill_ambiguous::Nullable{DNA}

    function FASTQReader(input::BufferedInputStream, quality_encoding, fill_ambiguous)
        return new(Ragel.State(fastqparser_start, input),
                   BufferedOutputStream(), BufferedOutputStream(),
                   StringField(), StringField(), 0, quality_encoding, fill_ambiguous)
    end
end

function Bio.IO.stream(reader::FASTQReader)
    return reader.state.stream
end

function Base.eltype{S}(::Type{FASTQReader{S}})
    return FASTQSeqRecord{S}
end

function FASTQReader(input::IO; quality_encoding=:sanger, fill_ambiguous=nothing)
    return FASTQReader{DNASequence}(BufferedInputStream(input), sym2qualenc(quality_encoding), fill_ambiguous)
end

function (::Type{FASTQReader{S}}){S<:Sequence}(input::IO;
                                               quality_encoding=:sanger,
                                               fill_ambiguous=nothing)
    return FASTQReader{S}(BufferedInputStream(input), sym2qualenc(quality_encoding), fill_ambiguous)
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
    elseif name == :illumina18
        return ILLUMINA18_QUAL_ENCODING
    end
    throw(ArgumentError("quality encoding ':$(name)' is not supported"))
end
