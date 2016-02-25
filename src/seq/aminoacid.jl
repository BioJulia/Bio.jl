# Aminoacids
# ==========

"Type representing AminoAcids"
bitstype 8 AminoAcid


# Conversion from/to integers
# ---------------------------

Base.convert(::Type{AminoAcid}, aa::UInt8) = reinterpret(AminoAcid, aa)
Base.convert(::Type{UInt8}, aa::AminoAcid) = reinterpret(UInt8, aa)
Base.convert{T<:Number}(::Type{T}, aa::AminoAcid) = convert(T, UInt8(aa))
Base.convert{T<:Number}(::Type{AminoAcid}, aa::T) = convert(AminoAcid, UInt8(aa))


# Arithmetic and Order
# --------------------

# These methods are necessary when deriving some algorithims
# like iteration, sort, comparison, and so on.
Base.(:-)(x::AminoAcid, y::AminoAcid) = Int(x) - Int(y)
Base.(:-)(x::AminoAcid, y::Integer) = reinterpret(AminoAcid, UInt8(x) - UInt8(y))
Base.(:+)(x::AminoAcid, y::Integer) = reinterpret(AminoAcid, UInt8(x) + UInt8(y))
Base.(:+)(x::Integer, y::AminoAcid) = y + x
Base.isless(x::AminoAcid, y::AminoAcid) = isless(UInt8(x), UInt8(y))


# Amino acid encoding definition
# ------------------------------

# This set of amino acids is defined by IUPAC-IUB Joint Commission on Biochemical Nomenclature.
# Reference: http://www.insdc.org/documents/feature_table.html#7.4.3

"Alanine"
const AA_A = convert(AminoAcid, 0x00)

"Arginine"
const AA_R = convert(AminoAcid, 0x01)

"Asparagine"
const AA_N = convert(AminoAcid, 0x02)

"Aspartic Acid"
const AA_D = convert(AminoAcid, 0x03)

"Cysteine"
const AA_C = convert(AminoAcid, 0x04)

"Glutamine"
const AA_Q = convert(AminoAcid, 0x05)

"Glutamic Acid"
const AA_E = convert(AminoAcid, 0x06)

"Glycine"
const AA_G = convert(AminoAcid, 0x07)

"Histidine"
const AA_H = convert(AminoAcid, 0x08)

"Isoleucine"
const AA_I = convert(AminoAcid, 0x09)

"Leucine"
const AA_L = convert(AminoAcid, 0x0a)

"Lysine"
const AA_K = convert(AminoAcid, 0x0b)

"Methionine"
const AA_M = convert(AminoAcid, 0x0c)

"Phenylalanine"
const AA_F = convert(AminoAcid, 0x0d)

"Proline"
const AA_P = convert(AminoAcid, 0x0e)

"Serine"
const AA_S = convert(AminoAcid, 0x0f)

"Threonine"
const AA_T = convert(AminoAcid, 0x10)

"Tryptophan"
const AA_W = convert(AminoAcid, 0x11)

"Tyrosine"
const AA_Y = convert(AminoAcid, 0x12)

"Valine"
const AA_V = convert(AminoAcid, 0x13)

"Aspartic Acid or Asparagine"  # ambiguous
const AA_B = convert(AminoAcid, 0x14)

"Leucine or Isoleucine"  # ambiguous
const AA_J = convert(AminoAcid, 0x15)

"Glutamine or Glutamic Acid"  # ambiguous
const AA_Z = convert(AminoAcid, 0x16)

"Unspecified or Unknown Amino Acid"  # ambiguous
const AA_X = convert(AminoAcid, 0x17)

"Pyrrolysine"  # non-standard
const AA_O = convert(AminoAcid, 0x18)

"Selenocysteine"  # non-standard
const AA_U = convert(AminoAcid, 0x19)

"Invalid Amino Acid"
const AA_INVALID = convert(AminoAcid, 0x1a) # Used during conversion from strings

isvalid(aa::AminoAcid) = aa ≤ AA_U
alphabet(::Type{AminoAcid}) = AA_A:AA_U

Base.isless(x::AminoAcid, y::AminoAcid) = isless(UInt8(x), UInt8(y))
Base.(:-)(x::AminoAcid, y::AminoAcid) = Int(UInt8(x)) - Int(UInt8(y))
Base.(:-)(x::AminoAcid, y::Integer) = reinterpret(AminoAcid, Int8(UInt8(x)) - Int8(UInt8(y)))
Base.(:+)(x::AminoAcid, y::Integer) = reinterpret(AminoAcid, Int8(UInt8(x)) + Int8(y))
Base.(:+)(x::Integer, y::AminoAcid) = y + x


# Conversion from/to Char
# -----------------------

# lookup table for characters in 'A':'z'
const char_to_aa = [
    AA_A,       AA_B,       AA_C,       AA_D,       AA_E,       AA_F,
    AA_G,       AA_H,       AA_I,       AA_J,       AA_K,       AA_L,
    AA_M,       AA_N,       AA_O,       AA_P,       AA_Q,       AA_R,
    AA_S,       AA_T,       AA_U,       AA_V,       AA_W,       AA_X,
    AA_Y,       AA_Z,       AA_INVALID, AA_INVALID, AA_INVALID, AA_INVALID,
    AA_INVALID, AA_INVALID, AA_A,       AA_B,       AA_C,       AA_D,
    AA_E,       AA_F,       AA_G,       AA_H,       AA_I,       AA_J,
    AA_K,       AA_L,       AA_M,       AA_N,       AA_O,       AA_P,
    AA_Q,       AA_R,       AA_S,       AA_T,       AA_U,       AA_V,
    AA_W,       AA_X,       AA_Y,       AA_Z]

function Base.convert(::Type{AminoAcid}, c::Char)
    @inbounds aa = 'A' <= c <= 'y' ? char_to_aa[c - 'A' + 1] : AA_INVALID
    @assert aa != AA_INVALID error("$(c) is not a valid amino acid")
    return aa
end

const aa_to_char = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
    'B', 'J', 'Z', 'X',
    'O', 'U', ]
Base.convert(::Type{Char}, aa::AminoAcid) = aa_to_char[convert(UInt8, aa) + 1]


# Basic functions
# ---------------

function Base.show(io::IO, aa::AminoAcid)
    if aa == AA_INVALID
        write(io, "Invalid Amino Acid")
    else
        write(io, convert(Char, aa))
    end
end

# lookup table of 20 standard amino acids
const threeletter_to_aa = Dict(
    "ALA" => AA_A, "ARG" => AA_R, "ASN" => AA_N, "ASP" => AA_D, "CYS" => AA_C,
    "GLN" => AA_Q, "GLU" => AA_E, "GLY" => AA_G, "HIS" => AA_H, "ILE" => AA_I,
    "LEU" => AA_L, "LYS" => AA_K, "MET" => AA_M, "PHE" => AA_F, "PRO" => AA_P,
    "SER" => AA_S, "THR" => AA_T, "TRP" => AA_W, "TYR" => AA_Y, "VAL" => AA_V,
    "ASX" => AA_B, "XLE" => AA_J, "GLX" => AA_Z, "XAA" => AA_X,
    "PYL" => AA_O, "SEC" => AA_U,
)

function Base.parse(::Type{AminoAcid}, s::AbstractString)
    s′ = strip(s)
    if length(s′) == 1
        return convert(AminoAcid, s′[1])
    end
    try
        return threeletter_to_aa[uppercase(s′)]
    catch ex
        if isa(ex, KeyError)
            error("invalid amino acid string: \"$s\" ")
        end
        rethrow()
    end
end

# TODO: tryparse
