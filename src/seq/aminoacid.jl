
bitstype 8 AminoAcid


# Conversion between integers and amino acids
function convert(::Type{AminoAcid}, aa::Uint8)
    return box(AminoAcid, unbox(Uint8, aa))
end


function convert(::Type{Uint8}, aa::AminoAcid)
    return box(Uint8, unbox(AminoAcid, aa))
end


function convert{T <: Unsigned}(::Type{T}, aa::AminoAcid)
    return box(T, Base.zext_int(T, unbox(AminoAcid, aa)))
end


function convert{T <: Unsigned}(::Type{AminoAcid}, aa::T)
    return convert(AminoAcid, convert(Uint8, aa))
end


# Amino acid encoding definition
const AA_A = convert(AminoAcid, 0x00)
const AA_R = convert(AminoAcid, 0x01)
const AA_N = convert(AminoAcid, 0x02)
const AA_D = convert(AminoAcid, 0x03)
const AA_C = convert(AminoAcid, 0x04)
const AA_Q = convert(AminoAcid, 0x05)
const AA_E = convert(AminoAcid, 0x06)
const AA_G = convert(AminoAcid, 0x07)
const AA_H = convert(AminoAcid, 0x08)
const AA_I = convert(AminoAcid, 0x09)
const AA_L = convert(AminoAcid, 0x0a)
const AA_K = convert(AminoAcid, 0x0b)
const AA_M = convert(AminoAcid, 0x0c)
const AA_F = convert(AminoAcid, 0x0d)
const AA_P = convert(AminoAcid, 0x0e)
const AA_S = convert(AminoAcid, 0x0f)
const AA_T = convert(AminoAcid, 0x10)
const AA_W = convert(AminoAcid, 0x11)
const AA_Y = convert(AminoAcid, 0x12)
const AA_V = convert(AminoAcid, 0x13)
const AA_X = convert(AminoAcid, 0x14)


# Use during conversion from strings
const AA_INVALID = convert(AminoAcid, 0x15)

const aa_to_char = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X' ]

function convert(::Type{Char}, aa::AminoAcid)
    return aa_to_char[convert(Uint8, aa) + 1]
end


# Conversion from Char

# lookup table for characters in 'A':'y'
const char_to_aa = [
    AA_A,       AA_INVALID, AA_C,       AA_D,       AA_E,       AA_F,
    AA_G,       AA_H,       AA_I,       AA_INVALID, AA_K,       AA_L,
    AA_M,       AA_N,       AA_INVALID, AA_P,       AA_Q,       AA_R,
    AA_S,       AA_T,       AA_INVALID, AA_V,       AA_W,       AA_X,
    AA_Y,       AA_INVALID, AA_INVALID, AA_INVALID, AA_INVALID, AA_INVALID,
    AA_INVALID, AA_INVALID, AA_A,       AA_INVALID, AA_C,       AA_D,
    AA_E,       AA_F,       AA_G,       AA_H,       AA_I,       AA_INVALID,
    AA_K,       AA_L,       AA_M,       AA_N,       AA_INVALID, AA_P,
    AA_Q,       AA_R,       AA_S,       AA_T,       AA_INVALID, AA_V,
    AA_W,       AA_X,       AA_Y ]


function convert(::Type{AminoAcid}, c::Char)
    @inbounds aa = 'A' <= c <= 'y' ? char_to_aa[c - 'A' + 1] : AA_INVALID
    if aa == AA_INVALID
        error("$(c) is not a valid amino acid")
    end

    return aa
end


function show(io::IO, aa::AminoAcid)
    write(io, convert(Char, aa))
end


# A general purpose amino acid representation.
#
# Amino acid are simple byte arrays using the encoding defined above. Any
# byte outside the range 0x00:0x14 is considered invalid and must result in an
# error.
#
# Like NucleotideSequence, amino acid sequences are immutable by convention.
#
type AminoAcidSequence
    data::Vector{AminoAcid}

    # interval within `data` defining the (sub)sequence
    part::UnitRange{Int}

    # Construct from raw components
    function AminoAcidSequence(data::Vector{AminoAcid}, part::UnitRange)
        return new(data, part)
    end

    # Construct a subsequence of another aa sequence
    function AminoAcidSequence(other::AminoAcidSequence, part::UnitRange)
        start = other.part.start + part.start - 1
        stop = start + length(part) - 1
        if start < other.part.start || stop > other.part.stop
            error("Invalid subsequence range")
        end
        return new(other.data, part)
    end

    # Construct of a subsequence from another amino acid sequence
    function AminoAcidSequence(seq::String)
        len = length(seq)
        data = Array(AminoAcid, len)
        for (i, c) in enumerate(seq)
            data[i] = convert(AminoAcid, c)
        end

        return new(data, 1:len)
    end
end


# Convert from/to String
function convert(::Type{AminoAcidSequence}, seq::String)
    return AminoAcidSequence(seq)
end


function convert(::Type{String}, seq::AminoAcidSequence)
    return convert(String, [convert(Char, x) for x in seq])
end


# String decorator syntax to enable building sequence literals like:
#     aa"ACDEFMN"
macro aa_str(seq, flags...)
    return AminoAcidSequence(seq)
end


function length(seq::AminoAcidSequence)
    return length(seq.part)
end


function endof(seq::AminoAcidSequence)
    return length(seq)
end


function show(io::IO, seq::AminoAcidSequence)
    const maxcount = 60
    write(io, "aa\"")
    len = length(seq)
    if len > maxcount
        for aa in seq[1:div(maxcount, 2) - 1]
            write(io, convert(Char, aa))
        end
        write("…")
        for aa in seq[(end - (div(maxcount, 2) - 1)):end]
            write(io, convert(Char, aa))
        end
    else
        for aa in seq
            write(io, convert(Char, aa))
        end
    end

    write(io, "\"  # ", string(len), "aa sequence")
end


function getindex(seq::NucleotideSequence, i::Integer)
    if i > length(seq) || i < 1
        error(BoundsError())
    end
    i += seq.part.start - 1
    return seq.data[i]
end


# Construct a subsequence
function getindex(seq::AminoAcidSequence, r::UnitRange)
    return AminoAcidSequence(seq, r)
end


# Replace a NucleotideSequence's data with a copy, copying only what's needed.
function orphan!(seq::AminoAcidSequence, reorphan=true)
    seq.data = seq.data[seq.part]
    seq.part = 1:length(seq.part)
    return seq
end


function copy(seq::AminoAcidSequence)
    return orphan!(AminoAcidSequence(seq.data, seq.part))
end


# Iterating through amino acids
function start(seq::AminoAcidSequence)
    return seq.part.start
end


function next(seq::AminoAcidSequence, i)
    aa = seq.data[i]
    return (aa, i + 1)
end


function done(seq::AminoAcidSequence, i)
    return i > length(seq)
end


# A genetic code is a table mapping RNA 3-mers (i.e. RNAKmer{3}) to AminoAcids.
immutable GeneticCode <: Associative{RNAKmer{3}, AminoAcid}
    tbl::Vector{AminoAcid}
end

function start(code::GeneticCode)
    return uint64(0)
end


function next(code::GeneticCode, x::Uint64)
    return ((convert(RNAKmer{3}, x), code[x]), (x + 1))
end


function done(code::GeneticCode, x::Uint64)
    return x > uint64(0b111111)
end


# All of these taken from:
# http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1
const standard_genetic_code = GeneticCode([
    AA_K, AA_Q, AA_E, AA_INVALID, AA_T, AA_P, AA_A, AA_S, AA_R, AA_R, AA_G,
    AA_INVALID, AA_I, AA_L, AA_V, AA_L, AA_N, AA_H, AA_D, AA_Y, AA_T, AA_P,
    AA_A, AA_S, AA_S, AA_R, AA_G, AA_C, AA_I, AA_L, AA_V, AA_F, AA_K, AA_Q,
    AA_E, AA_INVALID, AA_T, AA_P, AA_A, AA_S, AA_R, AA_R, AA_G, AA_W, AA_M,
    AA_L, AA_V, AA_L, AA_N, AA_H, AA_D, AA_Y, AA_T, AA_P, AA_A, AA_S, AA_S,
    AA_R, AA_G, AA_C, AA_I, AA_L, AA_V, AA_F ])

# TODO: more codes

function getindex(code::GeneticCode, idx::RNAKmer{3})
    return code.tbl[convert(Uint64, idx) + 1]
end


# Convert an RNASequence to an AminoAcidSequence
function translate(seq::RNASequence, code::GeneticCode=standard_genetic_code)
    d, r = divrem(length(seq), 3)
    if r != 0
        error("RNASequence length is not divisible by three. Cannot translate.")
    end

    if hasn(seq)
        error("Cannot translate an RNASequence with Ns")
    end

    aaseq = Array(AminoAcid, d)
    for (i, codon) in enumerate(eachkmer(seq, 3, 3))
        aa = code[codon]
        if aa == AA_INVALID
            error("Cannot translate stop codons.")
        end
        aaseq[i] = code[codon]
    end

    return AminoAcidSequence(aaseq, 1:d)
end




