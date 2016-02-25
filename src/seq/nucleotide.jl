# Nucleotides
# ===========


# Single nucleotides are represented in bytes using just the two low-order bits
abstract Nucleotide
bitstype 8 DNANucleotide <: Nucleotide
bitstype 8 RNANucleotide <: Nucleotide


# Conversion from/to integers
# ---------------------------

Base.convert(::Type{DNANucleotide}, nt::UInt8) = reinterpret(DNANucleotide, nt)
Base.convert(::Type{RNANucleotide}, nt::UInt8) = reinterpret(RNANucleotide, nt)
Base.convert(::Type{UInt8}, nt::DNANucleotide) = reinterpret(UInt8, nt)
Base.convert(::Type{UInt8}, nt::RNANucleotide) = reinterpret(UInt8, nt)
Base.convert{T<:Number,S<:Nucleotide}(::Type{T}, nt::S) = convert(T, UInt8(nt))
Base.convert{T<:Number,S<:Nucleotide}(::Type{S}, nt::T) = convert(S, UInt8(nt))


# Arithmetic and Order
# --------------------

# These methods are necessary when deriving some algorithims
# like iteration, sort, comparison, and so on.
Base.(:-){N<:Nucleotide}(x::N, y::N) = Int(x) - Int(y)
Base.(:-){N<:Nucleotide}(x::N, y::Integer) = reinterpret(N, UInt8(x) - UInt8(y))
Base.(:+){N<:Nucleotide}(x::N, y::Integer) = reinterpret(N, UInt8(x) + UInt8(y))
Base.(:+){N<:Nucleotide}(x::Integer, y::N) = y + x
Base.isless{N<:Nucleotide}(x::N, y::N) = isless(UInt8(x), UInt8(y))


# Nucleotide encoding definition
# ------------------------------

# DNA Nucleotides

"DNA Adenine"
const DNA_A = convert(DNANucleotide, 0b000)

"DNA Cytosine"
const DNA_C = convert(DNANucleotide, 0b001)

"DNA Guanine"
const DNA_G = convert(DNANucleotide, 0b010)

"DNA Thymine"
const DNA_T = convert(DNANucleotide, 0b011)

"DNA Any Nucleotide"
const DNA_N = convert(DNANucleotide, 0b100)

"DNA Invalid Nucleotide"
const DNA_INVALID = convert(DNANucleotide, 0b1000) # Indicates invalid DNA when converting string

"Returns Any DNA Nucleotide (DNA_N)"
nnucleotide(::Type{DNANucleotide}) = DNA_N
invalid_nucleotide(::Type{DNANucleotide}) = DNA_INVALID

isvalid(nt::DNANucleotide) = nt ≤ DNA_N
isambiguous(nt::DNANucleotide) = nt > DNA_T
alphabet(::Type{DNANucleotide}) = DNA_A:DNA_N

# RNA Nucleotides

"RNA Adenine"
const RNA_A = convert(RNANucleotide, 0b000)

"RNA Cytosine"
const RNA_C = convert(RNANucleotide, 0b001)

"RNA Guanine"
const RNA_G = convert(RNANucleotide, 0b010)

"RNA Uracil"
const RNA_U = convert(RNANucleotide, 0b011)

"Any RNA Nucleotide"
const RNA_N = convert(RNANucleotide, 0b100)

"Invalid RNA Nucleotide"
const RNA_INVALID = convert(RNANucleotide, 0b1000) # Indicates invalid RNA when converting string

"Returns Any RNA Nucleotide (RNA_N)"
nnucleotide(::Type{RNANucleotide}) = RNA_N
invalid_nucleotide(::Type{RNANucleotide}) = RNA_INVALID

isvalid(nt::RNANucleotide) = nt ≤ RNA_N
isambiguous(nt::RNANucleotide) = nt > RNA_U
alphabet(::Type{RNANucleotide}) = RNA_A:RNA_N


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

@inline function Base.convert(::Type{DNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 't' ? char_to_dna[c - 'A' + 1] : DNA_INVALID
    return nt
end

@inline function unsafe_ascii_byte_to_nucleotide(T::Type{DNANucleotide}, c::UInt8)
    @inbounds nt = char_to_dna[c - 0x41 + 1]
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

@inline function Base.convert(::Type{RNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 'u' ? char_to_rna[c - 'A' + 1] : RNA_INVALID
    return nt
end

@inline function unsafe_ascii_byte_to_nucleotide(T::Type{RNANucleotide}, c::UInt8)
    @inbounds nt = char_to_rna[c - 0x41 + 1]
    return nt
end


# Conversion to Char
# ------------------

const dna_to_char = ['A', 'C', 'G', 'T', 'N']

Base.convert(::Type{Char}, nt::DNANucleotide) = dna_to_char[convert(UInt8, nt) + 1]

const rna_to_char = ['A', 'C', 'G', 'U', 'N']

Base.convert(::Type{Char}, nt::RNANucleotide) = rna_to_char[convert(UInt8, nt) + 1]


# Basic functions
# ---------------

function Base.show(io::IO, nt::DNANucleotide)
    if !isvalid(nt)
        write(io, "Invalid DNA Nucleotide")
    else
        write(io, convert(Char, nt))
    end
end

function Base.show(io::IO, nt::RNANucleotide)
    if !isvalid(nt)
        write(io, "Invalid RNA Nucleotide")
    else
        write(io, convert(Char, nt))
    end
end
