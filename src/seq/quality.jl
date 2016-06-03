# Base Quality
# ============
#
# A representation of positions-specific integer quality scores, as in FASTQ.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"A `QualityEncoding` value holds a set of encoding compatible with a string"
bitstype 16 QualityEncoding

Base.convert(::Type{QualityEncoding}, nt::UInt16) = reinterpret(QualityEncoding, nt)
Base.convert(::Type{UInt16}, nt::QualityEncoding) = reinterpret(UInt16, nt)

@compat function Base.:|(a::QualityEncoding, b::QualityEncoding)
    return convert(QualityEncoding, convert(UInt16, a) | convert(UInt16, b))
end

@compat function Base.:&(a::QualityEncoding, b::QualityEncoding)
    return convert(QualityEncoding, convert(UInt16, a) & convert(UInt16, b))
end

"`QualityEncoding` indicating no compatible encodings."
const EMPTY_QUAL_ENCODING      = convert(QualityEncoding, UInt16(0))

"Sanger (Phred+33) quality score encoding"
const SANGER_QUAL_ENCODING     = convert(QualityEncoding, UInt16(0b00001))

"Solexa (Solexa+64) quality score encoding"
const SOLEXA_QUAL_ENCODING     = convert(QualityEncoding, UInt16(0b00010))

"Illumina 1.3 (Phred+64) quality score encoding"
const ILLUMINA13_QUAL_ENCODING = convert(QualityEncoding, UInt16(0b00100))

"Illumina 1.5 (Phred+64) quality score encoding"
const ILLUMINA15_QUAL_ENCODING = convert(QualityEncoding, UInt16(0b01000))

"Illumina 1.8 (Phred+33) quality score encoding"
const ILLUMINA18_QUAL_ENCODING = convert(QualityEncoding, UInt16(0b10000))

"`QualityEncoding` indicating all known encodings are compatible."
const ALL_QUAL_ENCODINGS =
    SANGER_QUAL_ENCODING | SOLEXA_QUAL_ENCODING | ILLUMINA13_QUAL_ENCODING |
    ILLUMINA15_QUAL_ENCODING | ILLUMINA18_QUAL_ENCODING

# Index into this with `trailing_zeros(encoding) + 1`
const qual_encoding_ranges = [
    (UInt8('!'), UInt8('~'), Int8(0) ),  # SANGER
    (UInt8(';'), UInt8('~'), Int8(-5)),  # SOLEXA
    (UInt8('@'), UInt8('~'), Int8(0) ),  # ILLUMINA13
    (UInt8('B'), UInt8('~'), Int8(2) ),  # ILLUMINA15
    (UInt8('!'), UInt8('~'), Int8(0) )   # ILLUMINA18
]

# Build an encoding lookup table
const compatible_qual_encoding = fill(EMPTY_QUAL_ENCODING, length('!':'~'))

for (encoding_num, (first, last, start_score)) in enumerate(qual_encoding_ranges)
    encoding = convert(QualityEncoding, UInt16(1) << (encoding_num - 1))
    for c in first:last
        compatible_qual_encoding[c - UInt8('!') + 1] |= encoding
    end
end

# When a quality string has multiple compatible encodings, we choose the first
# compatible alphabet in this list.
const preferred_quality_encodings = [
    ILLUMINA15_QUAL_ENCODING,
    ILLUMINA13_QUAL_ENCODING,
    SOLEXA_QUAL_ENCODING,
    SANGER_QUAL_ENCODING,
    ILLUMINA18_QUAL_ENCODING,
]

"""
Infer the encoding of a FASTQ quality string.

# Arguments
  * `data`: sequence data in a string
  * `start`: first position to consider in data
  * `stop`: last position to consider in data
  * `encodings`: valid encodings
  * `default`: if there are multiple compatible encodings, default
    to this one if it's compatible.

# Returns
A pair of `QualityEncoding`s. The first given the chosen encoding, the second
giving the set of compatible encodings.
"""
function infer_quality_encoding(data::Vector{UInt8}, start, stop,
                                encodings::QualityEncoding=ALL_QUAL_ENCODINGS,
                                default::QualityEncoding=EMPTY_QUAL_ENCODING)
    encodings = ALL_QUAL_ENCODINGS
    @inbounds for i in start:stop
        c = data[i]
        if '!' <= c <= '~'
            encodings &= compatible_qual_encoding[c - UInt8('!') + 1]
        else
            error("Character $(convert(Char, c)) is not compatible with any known quality encoding.")
        end

        if count_ones(convert(UInt16, encodings)) <= 0
            break
        end
    end

    if count_ones(convert(UInt16, encodings)) == 0
        error("String is not compatible with any known sequence type.")
    elseif default != EMPTY_QUAL_ENCODING && (encodings & default) != EMPTY_ALPHABET
        return (default, encodings)
    else
        for encoding in preferred_quality_encodings
            if encoding & encodings != EMPTY_QUAL_ENCODING
                return (encoding, encodings)
            end
        end
    end
    return (default, encodings)
end

"""
Decode a quality string in place into integer Phred scores.

# Arguments:
  * `encoding::QualityEncoding`: how the quality scores are encoded.
  * `input::Vector{UInt8}`: character data to decode
  * `output:::Vector{UInt8}`: Phred score output vector, assumed to be
    of length `stop - start + 1`.
  * `start`: First position in `input` to decode.
  * `stop`: Last position with in `input` to decode.

# Returns
`output`
"""
function decode_quality_string!(encoding::QualityEncoding, input::Vector{UInt8},
                                output::Vector{Int8}, start=1,
                                stop=min(length(output), length(input)))
    if length(output) != stop - start + 1
        resize!(output, stop - start + 1)
    end

    encoding_num = trailing_zeros(convert(UInt16, encoding)) + 1
    first, last, startqual = qual_encoding_ranges[encoding_num]

    fill!(output, startqual - first)
    for i in 1:(stop - start + 1)
        @inbounds output[i] += input[start + i - 1]
    end

    return output
end

"""
Decode a quality string in place into integer Phred scores.

# Arguments:
  * `encoding::QualityEncoding`: how the quality scores are encoded.
  * `input::Vector{UInt8}`: character data to decode
  * `start`: First position in `input` to decode.
  * `stop`: Last position with in `input` to decode.

# Returns
A `Vector{UInt8}` of length `stop - start + 1` containing integer Phred scores.
"""
function decode_quality_string(encoding::QualityEncoding, input::Vector{UInt8},
                               start=1, stop=length(input))
    output = Array(Int8, stop - start + 1)
    return decode_quality_string!(encoding, input, output, start, stop)
end

"""
Encode a Phred quality score vector in place into a quality string.

# Arguments:
  * `encoding::QualityEncoding`: how the quality scores are encoded.
  * `input::Vector{UInt8}`: character data to decode
  * `output:::Vector{UInt8}`: Phred score output vector, assumed to be
    of length `stop - start + 1`.
  * `start`: First position in `input` to decode.
  * `stop`: Last position with in `input` to decode.

# Returns
`output`
"""
function encode_quality_string!(encoding::QualityEncoding, input::Vector{Int8},
                                output::Vector{UInt8}, start=1,
                                stop=min(length(output), length(input)))
    @inbounds begin
        if length(output) != stop - start + 1
            resize!(output, stop - start + 1)
        end

        encoding_num = trailing_zeros(convert(UInt16, encoding)) + 1
        first, last, startqual = qual_encoding_ranges[encoding_num]
        for (i, j) in enumerate(start:stop)
            c = input[j]
            output[i] =  (c + first) - startqual
        end
        return output
    end
end

"""
Encode a vector of integer Phred scores into a quality string

# Arguments:
  * `encoding::QualityEncoding`: how the quality scores are encoded.
  * `input::Vector{UInt8}`: character data to decode
  * `start`: First position in `input` to decode.
  * `stop`: Last position with in `input` to decode.

# Returns
A `Vector{UInt8}` of length `stop - start + 1` containing ASCII-encoded phred
scores
"""
function encode_quality_string(encoding::QualityEncoding, input::Vector{Int8},
                               start=1, stop=length(input))
    output = Array(UInt8, stop - start + 1)
    return encode_quality_string!(encoding, input, output, start, stop)
end
