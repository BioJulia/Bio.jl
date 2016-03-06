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

immutable EncodeError{A<:Alphabet,T} <: Exception
    val::T
end

Base.call{A,T}(::Type{EncodeError{A}}, val::T) = EncodeError{A,T}(val)

function Base.showerror{A}(io::IO, err::EncodeError{A})
    print(io, "cannot encode ", err.val, " in ", A)
end

"""
Decode biological characters from binary representation.
"""
function decode end

immutable DecodeError{A<:Alphabet,T} <: Exception
    val::T
end

Base.call{A,T}(::Type{DecodeError{A}}, val::T) = DecodeError{A,T}(val)

function Base.showerror{A}(io::IO, err::DecodeError{A})
    print(io, "cannot decode ", err.val, " in ", A)
end

for (A, T, ub) in [
        (DNAAlphabet{2},    DNANucleotide, DNA_T  ),
        (DNAAlphabet{4},    DNANucleotide, DNA_Gap),
        (RNAAlphabet{2},    RNANucleotide, RNA_U  ),
        (RNAAlphabet{4},    RNANucleotide, RNA_Gap),
        (AminoAcidAlphabet, AminoAcid,     AA_Gap )]
    @eval begin
        @inline function encode(::Type{$A}, x::$T)
            if x > $ub
                throw(EncodeError{$A}(x))
            end
            return reinterpret(UInt8, x)
        end
        @inline function decode(::Type{$A}, x::UInt8)
            if x > $(reinterpret(UInt8, ub))
                throw(DecodeError{$A}(x))
            end
            return reinterpret($T, x)
        end
    end
end
