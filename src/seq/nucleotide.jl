# Nucleotides
# ===========
#
# DNA and RNA nucleotide types.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Valid nucleotides are represented in bytes using just the four low-order bits
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
@compat begin
    Base.:-{N<:Nucleotide}(x::N, y::N) = Int(x) - Int(y)
    Base.:-{N<:Nucleotide}(x::N, y::Integer) = reinterpret(N, UInt8(x) - UInt8(y))
    Base.:+{N<:Nucleotide}(x::N, y::Integer) = reinterpret(N, UInt8(x) + UInt8(y))
    Base.:+{N<:Nucleotide}(x::Integer, y::N) = y + x
end
Base.isless{N<:Nucleotide}(x::N, y::N) = isless(UInt8(x), UInt8(y))


# Nucleotide encoding definition
# ------------------------------

# DNA Nucleotides

"DNA Invalid Nucleotide"
const DNA_INVALID = convert(DNANucleotide, 0b10000) # Indicates invalid DNA when converting string

# lookup table for characters
const char_to_dna = [DNA_INVALID for _ in 0x00:0x7f]
const dna_to_char = Vector{Char}(16)

# compatibility bits
const compatbits = zeros(UInt8, 16)

# derived from "The DDBJ/ENA/GenBank Feature Table Definition"
# §7.4.1 Nucleotide base code (IUPAC)
# http://www.insdc.org/documents/feature_table.html#7.4.1
for (code, (nt, doc, compat)) in enumerate([
        ('A', "DNA Adenine",                               0b0001),
        ('C', "DNA Cytosine",                              0b0010),
        ('G', "DNA Guanine",                               0b0100),
        ('T', "DNA Thymine",                               0b1000),
        ('M', "DNA Adenine or Cytosine",                   0b0011),
        ('R', "DNA Adenine or Guanine",                    0b0101),
        ('W', "DNA Adenine or Thymine",                    0b1001),
        ('S', "DNA Cytosine or Guanine",                   0b0110),
        ('Y', "DNA Cytosine or Thymine",                   0b1010),
        ('K', "DNA Guanine or Thymine",                    0b1100),
        ('V', "DNA Adenine, Cytosine or Guanine",          0b0111),
        ('H', "DNA Adenine, Cytosine or Thymine",          0b1011),
        ('D', "DNA Adenine, Guanine or Thymine",           0b1101),
        ('B', "DNA Cytosine, Guanine or Thymine",          0b1110),
        ('N', "DNA Adenine, Cytosine, Guanine or Thymine", 0b1111)])
    var = Symbol("DNA_", nt)
    @eval begin
        @doc $doc const $var = convert(DNANucleotide, $(UInt8(code - 1)))
        char_to_dna[$(Int(nt + 1))] = char_to_dna[$(Int(lowercase(nt) + 1))] = $var
        dna_to_char[$(code)] = $nt
        compatbits[$(code)] = $compat
    end
end

"DNA Gap"
const DNA_Gap = convert(DNANucleotide, 0b1111)
char_to_dna[Int('-') + 1] = DNA_Gap
dna_to_char[0b1111 + 1] = '-'
compatbits[0b1111 + 1] = 0b0000

"Returns Any DNA Nucleotide (DNA_N)"
nnucleotide(::Type{DNANucleotide}) = DNA_N
invalid_nucleotide(::Type{DNANucleotide}) = DNA_INVALID

Base.isvalid(::Type{DNANucleotide}, x::Integer) = 0 ≤ x < 16
Base.isvalid(nt::DNANucleotide) = nt ≤ DNA_Gap
isambiguous(nt::DNANucleotide) = nt > DNA_T
alphabet(::Type{DNANucleotide}) = DNA_A:DNA_Gap
gap(::Type{DNANucleotide}) = DNA_Gap

# check if two nucleotides are compatible (no contradiction) with each other
function iscompatible(x::Nucleotide, y::Nucleotide)
    return (compatbits[Int8(x)+1] & compatbits[Int8(y)+1]) != 0
end

# RNA Nucleotides

"Invalid RNA Nucleotide"
const RNA_INVALID = convert(RNANucleotide, 0b10000) # Indicates invalid RNA when converting string

# lookup table for characters
const char_to_rna = [RNA_INVALID for _ in 0x00:0x7f]
const rna_to_char = Vector{Char}(16)

for (code, (nt, doc)) in enumerate([
        ('A', "RNA Adenine"                             ),
        ('C', "RNA Cytosine"                            ),
        ('G', "RNA Guanine"                             ),
        ('U', "RNA Uracil"                              ),
        ('M', "RNA Adenine or Cytosine"                 ),
        ('R', "RNA Adenine or Guanine"                  ),
        ('W', "RNA Adenine or Uracil"                   ),
        ('S', "RNA Cytosine or Guanine"                 ),
        ('Y', "RNA Cytosine or Uracil"                  ),
        ('K', "RNA Guanine or Uracil"                   ),
        ('V', "RNA Adenine, Cytosine or Guanine"        ),
        ('H', "RNA Adenine, Cytosine or Uracil"         ),
        ('D', "RNA Adenine, Guanine or Uracil"          ),
        ('B', "RNA Cytosine, Guanine or Uracil"         ),
        ('N', "RNA Adenine, Cytosine, Guanine or Uracil")])
    var = Symbol("RNA_", nt)
    @eval begin
        @doc $doc const $var = convert(RNANucleotide, $(UInt8(code - 1)))
        char_to_rna[$(Int(nt + 1))] = char_to_rna[$(Int(lowercase(nt) + 1))] = $var
        rna_to_char[$(code)] = $nt
    end
end

"RNA Gap"
const RNA_Gap = convert(RNANucleotide, 0b1111)
char_to_rna[Int('-') + 1] = RNA_Gap
rna_to_char[0b1111 + 1] = '-'

"Returns Any RNA Nucleotide (RNA_N)"
nnucleotide(::Type{RNANucleotide}) = RNA_N
invalid_nucleotide(::Type{RNANucleotide}) = RNA_INVALID

Base.isvalid(::Type{RNANucleotide}, x::Integer) = 0 ≤ x < 16
Base.isvalid(nt::RNANucleotide) = nt ≤ RNA_Gap
isambiguous(nt::RNANucleotide) = nt > RNA_U
alphabet(::Type{RNANucleotide}) = RNA_A:RNA_Gap
gap(::Type{RNANucleotide}) = RNA_Gap


# Conversion from Char
# --------------------

function Base.convert(::Type{DNANucleotide}, c::Char)
    @inbounds return c <= '\x7f' ? char_to_dna[Int(c) + 1] : DNA_INVALID
end

function unsafe_ascii_byte_to_nucleotide(T::Type{DNANucleotide}, c::UInt8)
    @inbounds return char_to_dna[c + 1]
end

function Base.convert(::Type{RNANucleotide}, c::Char)
    @inbounds return c <= '\x7f' ? char_to_rna[Int(c) + 1] : RNA_INVALID
end

function unsafe_ascii_byte_to_nucleotide(T::Type{RNANucleotide}, c::UInt8)
    @inbounds return char_to_rna[c + 1]
end


# Conversion to Char
# ------------------

Base.convert(::Type{Char}, nt::DNANucleotide) = dna_to_char[convert(UInt8, nt) + 1]
Base.convert(::Type{Char}, nt::RNANucleotide) = rna_to_char[convert(UInt8, nt) + 1]


# Basic functions
# ---------------

function Base.show(io::IO, nt::DNANucleotide)
    if !isvalid(nt)
        write(io, "Invalid DNA Nucleotide")
    else
        write(io, Char(nt))
    end
end

function Base.show(io::IO, nt::RNANucleotide)
    if !isvalid(nt)
        write(io, "Invalid RNA Nucleotide")
    else
        write(io, Char(nt))
    end
end
