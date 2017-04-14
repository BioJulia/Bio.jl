# FASTQ Base Quality
# ==================
#
# A representation of positions-specific integer quality scores, as in FASTQ.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable QualityEncoding
    # valid ASCII character range
    ascii::NTuple{2,UInt8}

    # valid base quality range
    qual::NTuple{2,Int8}

    function QualityEncoding(ascii::Range{Char}, qual::Range{Int})
        if length(ascii) != length(qual)
            throw(ArgumentError("the ranges of ASCII and base quality don't match"))
        end
        return new((first(ascii), last(ascii)), (first(qual), last(qual)))
    end
end

# Return the offset of converting a base quality to an ASCII code (quality + offset).
function ascii_offset(encoding::QualityEncoding)
    return encoding.ascii[1] - encoding.qual[1]
end

"Sanger (Phred+33) quality score encoding"
const SANGER_QUAL_ENCODING     = QualityEncoding('!':'~', 0:93)

"Solexa (Solexa+64) quality score encoding"
const SOLEXA_QUAL_ENCODING     = QualityEncoding(';':'~', -5:62)

"Illumina 1.3 (Phred+64) quality score encoding"
const ILLUMINA13_QUAL_ENCODING = QualityEncoding('@':'~', 0:62)

"Illumina 1.5 (Phred+64) quality score encoding"
const ILLUMINA15_QUAL_ENCODING = QualityEncoding('B':'~', 2:62)

"Illumina 1.8 (Phred+33) quality score encoding"
const ILLUMINA18_QUAL_ENCODING = QualityEncoding('!':'~', 0:93)


# Check quality encoding of `input[start:stop]`.
function check_quality_string(encoding, input, start, stop)
    ascii_lo, ascii_hi = encoding.ascii
    for i in start:stop
        @inbounds a = input[i]
        if !(ascii_lo ≤ a ≤ ascii_hi)
            error("base quality '$(Char(a))' is out of range")
        end
    end
end

# Decode an ASCII string `input[start:stop]` into base quality scores `output`.
function decode_quality_string!(encoding, input, output, start, stop)
    resize!(output, stop - start + 1)
    offset = ascii_offset(encoding)
    for i in 1:(stop - start + 1)
        @inbounds output[i] = input[start + i - 1] - offset
    end
    return output
end

# Encode base quality scores `input` into a ASCII string `output`.
function encode_quality_string!(encoding, input, output, start, stop)
    resize!(output, stop - start + 1)
    offset = ascii_offset(encoding)
    for i in 1:(stop - start + 1)
        @inbounds output[i] = input[start + i - 1] + offset
    end
    return output
end
