# Aminoacids
# ==========


bitstype 8 AminoAcid


# Conversion from/to integers
# ---------------------------

convert(::Type{AminoAcid}, aa::Uint8) = box(AminoAcid, unbox(Uint8, aa))
convert(::Type{Uint8}, aa::AminoAcid) = box(Uint8, unbox(AminoAcid, aa))
convert{T <: Unsigned}(::Type{T}, aa::AminoAcid) = box(T, Base.zext_int(T, unbox(AminoAcid, aa)))
convert{T <: Unsigned}(::Type{AminoAcid}, aa::T) = convert(AminoAcid, convert(Uint8, aa))


# Amino acid encoding definition
# ------------------------------

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
const AA_INVALID = convert(AminoAcid, 0x15) # Used during conversion from strings


# Conversion from/to Char
# -----------------------

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
    @assert aa != AA_INVALID error("$(c) is not a valid aminoacid")
    return aa
end

const aa_to_char = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X' ]

convert(::Type{Char}, aa::AminoAcid) = aa_to_char[convert(Uint8, aa) + 1]


# Basic functions
# ---------------

function show(io::IO, aa::AminoAcid)
    if aa == AA_INVALID
        write(io, "Invalid Amino Acid")
    else
        write(io, convert(Char, aa))
    end
end
