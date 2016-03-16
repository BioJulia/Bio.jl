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

"Invalid Amino Acid"
const AA_INVALID = convert(AminoAcid, 0x1c)  # Used during conversion from strings

# lookup table for characters
const char_to_aa = [AA_INVALID for _ in 0x00:0x7f]
const aa_to_char = Vector{Char}(0x1c)

# This set of amino acids is defined by IUPAC-IUB Joint Commission on Biochemical Nomenclature.
# Reference: http://www.insdc.org/documents/feature_table.html#7.4.3

for (aa, doc, code) in [
        ('A', "Alanine",                           0x00),
        ('R', "Arginine",                          0x01),
        ('N', "Asparagine",                        0x02),
        ('D', "Aspartic Acid",                     0x03),
        ('C', "Cysteine",                          0x04),
        ('Q', "Glutamine",                         0x05),
        ('E', "Glutamic Acid",                     0x06),
        ('G', "Glycine",                           0x07),
        ('H', "Histidine",                         0x08),
        ('I', "Isoleucine",                        0x09),
        ('L', "Leucine",                           0x0a),
        ('K', "Lysine",                            0x0b),
        ('M', "Methionine",                        0x0c),
        ('F', "Phenylalanine",                     0x0d),
        ('P', "Proline",                           0x0e),
        ('S', "Serine",                            0x0f),
        ('T', "Threonine",                         0x10),
        ('W', "Tryptophan",                        0x11),
        ('Y', "Tyrosine",                          0x12),
        ('V', "Valine",                            0x13),
        ('B', "Aspartic Acid or Asparagine",       0x14),  # ambiguous
        ('J', "Leucine or Isoleucine",             0x15),  # ambiguous
        ('Z', "Glutamine or Glutamic Acid",        0x16),  # ambiguous
        ('X', "Unspecified or Unknown Amino Acid", 0x17),  # ambiguous
        ('O', "Pyrrolysine",                       0x18),  # non-standard
        ('U', "Selenocysteine",                    0x19)]  # non-standard
    var = symbol("AA_", aa)
    @eval begin
        @doc $doc const $var = convert(AminoAcid, $code)
        char_to_aa[$(Int(aa)+1)] = char_to_aa[$(Int(lowercase(aa))+1)] = $var
        aa_to_char[$(code)+1] = $aa
    end
end

"Stop"
const AA_Stop = convert(AminoAcid, 0x1a)
char_to_aa[Int('*')+1] = AA_Stop
aa_to_char[0x1a+1] = '*'

"Amino Acid Gap"
const AA_Gap = convert(AminoAcid, 0x1b)
char_to_aa[Int('-') + 1] = AA_Gap
aa_to_char[0x1b+1] = '-'

Base.isvalid(::Type{AminoAcid}, x::Integer) = 0 ≤ x ≤ 0x1b
Base.isvalid(aa::AminoAcid) = aa ≤ AA_Gap
alphabet(::Type{AminoAcid}) = AA_A:AA_Gap


# Conversion from/to Char
# -----------------------

function Base.convert(::Type{AminoAcid}, c::Char)
    @inbounds aa = c <= '\x7f' ? char_to_aa[Int(c)+1] : AA_INVALID
    @assert aa != AA_INVALID error("$(c) is not a valid amino acid")
    return aa
end

Base.convert(::Type{Char}, aa::AminoAcid) = aa_to_char[convert(UInt8, aa) + 1]


# Basic functions
# ---------------

function Base.show(io::IO, aa::AminoAcid)
    if aa == AA_INVALID
        write(io, "Invalid Amino Acid")
    else
        write(io, Char(aa))
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
