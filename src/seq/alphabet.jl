# Alphabet
# ========
#
# Subtypes of `Alphabet` represent a domain of biological characters. For
# example, `DNAAlphabet{2}` has a domain of unambiguous nucleotides (i.e. A, C,
# G, and T). These types are used for parameterizing biological sequences and so
# on. A pair of encoder and decoder is associated with an alphabet, which maps
# values between binary and Julia-level representation.

"""
Alphabet of biological characters.
"""
abstract Alphabet

"""
DNA nucleotide alphabet.
"""
immutable DNAAlphabet{n} <: Alphabet end

"""
RNA nucleotide alphabet.
"""
immutable RNAAlphabet{n} <: Alphabet end

"""
Amino acid alphabet.
"""
immutable AminoAcidAlphabet <: Alphabet end

"""
The number of bits to represent the alphabet.
"""
function bitsof end

for n in (2, 4)
    @eval begin
        bitsof(::Type{DNAAlphabet{$n}}) = $n
        bitsof(::Type{RNAAlphabet{$n}}) = $n
    end
end
bitsof(::Type{AminoAcidAlphabet}) = 8

Base.eltype(::Type{DNAAlphabet}) = DNANucleotide
Base.eltype(::Type{RNAAlphabet}) = RNANucleotide
Base.eltype{n}(::Type{DNAAlphabet{n}}) = DNANucleotide
Base.eltype{n}(::Type{RNAAlphabet{n}}) = RNANucleotide
Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid

alphabet(::Type{DNAAlphabet{2}}) = DNA_A:DNA_T
alphabet(::Type{RNAAlphabet{2}}) = RNA_A:RNA_U
alphabet(::Type{DNAAlphabet{4}}) = alphabet(DNANucleotide)
alphabet(::Type{RNAAlphabet{4}}) = alphabet(RNANucleotide)
alphabet(::Type{AminoAcidAlphabet}) = alphabet(AminoAcid)


# Encoders & Decoders
# -------------------

"""
Encode biological characters to binary representation.
"""
function encode end

"""
Decode biological characters from binary representation.
"""
function decode end

for (A, N) in ((DNAAlphabet, DNANucleotide),
               (RNAAlphabet, RNANucleotide)), n in (2, 4)
    @eval begin
        encode(::Type{$A{$n}}, x::$N) = reinterpret(UInt8, x)
        decode(::Type{$A{$n}}, x::UInt8) = reinterpret($N, x)
    end
end

encode(::Type{AminoAcidAlphabet}, x::AminoAcid) = reinterpret(UInt8, x)
decode(::Type{AminoAcidAlphabet}, x::UInt8) = reinterpret(AminoAcid, x)
