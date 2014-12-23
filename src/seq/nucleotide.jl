# Nucleotides
# ===========


# Single nucleotides are represented in bytes using just the two low-order bits
abstract Nucleotide

bitstype 8 DNANucleotide <: Nucleotide
bitstype 8 RNANucleotide <: Nucleotide


# Conversion from/to integers
# ---------------------------

convert(::Type{DNANucleotide}, nt::Uint8) = box(DNANucleotide, unbox(Uint8, nt))
convert(::Type{Uint8}, nt::DNANucleotide) = box(Uint8, unbox(DNANucleotide, nt))

convert(::Type{RNANucleotide}, nt::Uint8) = box(RNANucleotide, unbox(Uint8, nt))
convert(::Type{Uint8}, nt::RNANucleotide) = box(Uint8, unbox(RNANucleotide, nt))

convert{T<:Unsigned, S<:Nucleotide}(::Type{T}, nt::S) = box(T, Base.zext_int(T, unbox(S, nt)))
convert{T<:Unsigned, S<:Nucleotide}(::Type{S}, nt::T) = convert(S, convert(Uint8, nt))


# Nucleotide encoding definition
# ------------------------------

const DNA_A = convert(DNANucleotide, 0b000)
const DNA_C = convert(DNANucleotide, 0b001)
const DNA_G = convert(DNANucleotide, 0b010)
const DNA_T = convert(DNANucleotide, 0b011)
const DNA_N = convert(DNANucleotide, 0b100)
const DNA_INVALID = convert(DNANucleotide, 0b1000) # Indicates invalid DNA when converting string

nnucleotide(::Type{DNANucleotide}) = DNA_N

const RNA_A = convert(RNANucleotide, 0b000)
const RNA_C = convert(RNANucleotide, 0b001)
const RNA_G = convert(RNANucleotide, 0b010)
const RNA_U = convert(RNANucleotide, 0b011)
const RNA_N = convert(RNANucleotide, 0b100)
const RNA_INVALID = convert(RNANucleotide, 0b1000) # Indicates invalid DNA when converting string

nnucleotide(::Type{RNANucleotide}) = RNA_N


# Conversion from Char
# --------------------

# lookup table for characters in 'A':'t'
const char_to_dna = [
     DNA_A,       DNA_INVALID, DNA_C,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_N,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_T,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_A,       DNA_INVALID, DNA_C,       DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_N,       DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_T ]


function convert(::Type{DNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 't' ? char_to_dna[c - 'A' + 1] : DNA_INVALID
    @assert nt != DNA_INVALID error("$(c) is not a valid DNA nucleotide")
    return nt
end

# lookup table for characters in 'A':'u'
const char_to_rna = [
    RNA_A,       RNA_INVALID, RNA_C,       RNA_INVALID, RNA_INVALID, RNA_INVALID,
    RNA_G,       RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_INVALID,
    RNA_INVALID, RNA_N,       RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_INVALID,
    RNA_INVALID, RNA_INVALID, RNA_U,       RNA_INVALID, RNA_INVALID, RNA_INVALID,
    RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_INVALID,
    RNA_INVALID, RNA_INVALID, RNA_A,       RNA_INVALID, RNA_C,       RNA_INVALID,
    RNA_INVALID, RNA_INVALID, RNA_G,       RNA_INVALID, RNA_INVALID, RNA_INVALID,
    RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_N,       RNA_INVALID, RNA_INVALID,
    RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_INVALID, RNA_U ]

function convert(::Type{RNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 'u' ? char_to_rna[c - 'A' + 1] : RNA_INVALID
    @assert nt != RNA_INVALID error("$(c) is not a valid RNA nucleotide")
    return nt
end


# Conversion to Char
# ------------------

const dna_to_char = ['A', 'C', 'G', 'T', 'N']
convert(::Type{Char}, nt::DNANucleotide) = dna_to_char[convert(Uint8, nt) + 1]

const rna_to_char = ['A', 'C', 'G', 'U', 'N']
convert(::Type{Char}, nt::RNANucleotide) = rna_to_char[convert(Uint8, nt) + 1]


# Basic functions
# ---------------

function show{T<:Nucleotide}(io::IO, nt::T)
    if nt == DNA_INVALID
        write(io, "Invalid DNA Nucleotide")
    elseif nt == RNA_INVALID
        write(io, "Invalid RNA Nucleotide")
    else
        write(io, convert(Char, nt))
    end
end
