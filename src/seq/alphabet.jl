# Alphabet
# ========
#
# Alphabet of biological symbols.
#
# Subtypes of `Alphabet` represent a domain of biological characters. For
# example, `DNAAlphabet{2}` has a domain of unambiguous nucleotides (i.e. A, C,
# G, and T). These types are used for parameterizing biological sequences and so
# on. A pair of encoder and decoder is associated with an alphabet, which maps
# values between binary and Julia-level representation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

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
General character alphabet.
"""
immutable CharAlphabet <: Alphabet end

"""
Void alphabet (internal use only).
"""
immutable VoidAlphabet <: Alphabet end

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
bitsof(::Type{CharAlphabet}) = 32
bitsof(::Type{VoidAlphabet}) = 0

Base.eltype(::Type{DNAAlphabet}) = DNANucleotide
Base.eltype(::Type{RNAAlphabet}) = RNANucleotide
Base.eltype{n}(::Type{DNAAlphabet{n}}) = DNANucleotide
Base.eltype{n}(::Type{RNAAlphabet{n}}) = RNANucleotide
Base.eltype(::Type{AminoAcidAlphabet}) = AminoAcid
Base.eltype(::Type{CharAlphabet}) = Char
Base.eltype(::Type{VoidAlphabet}) = Void

alphabet(::Type{DNAAlphabet{2}}) = DNA_A:DNA_T
alphabet(::Type{RNAAlphabet{2}}) = RNA_A:RNA_U
alphabet(::Type{DNAAlphabet{4}}) = alphabet(DNANucleotide)
alphabet(::Type{RNAAlphabet{4}}) = alphabet(RNANucleotide)
alphabet(::Type{AminoAcidAlphabet}) = alphabet(AminoAcid)
# TODO: this alphabet includes invalid Unicode scalar values
alphabet(::Type{CharAlphabet}) = typemin(Char):typemax(Char)
alphabet(::Type{VoidAlphabet}) = nothing


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

for (A, T, U, ub) in [
        (DNAAlphabet{2},    DNANucleotide, UInt8,  DNA_T     ),
        (DNAAlphabet{4},    DNANucleotide, UInt8,  DNA_Gap   ),
        (RNAAlphabet{2},    RNANucleotide, UInt8,  RNA_U     ),
        (RNAAlphabet{4},    RNANucleotide, UInt8,  RNA_Gap   ),
        (AminoAcidAlphabet, AminoAcid,     UInt8,  AA_Gap    ),
        (CharAlphabet,      Char,          UInt32, '\U10ffff')]
    @assert sizeof(T) == sizeof(U)
    @assert isa(ub, T)
    @eval begin
        @inline function encode(::Type{$A}, x::$T)
            if x > $ub
                throw(EncodeError{$A}(x))
            end
            return reinterpret($U, x)
        end
        @inline function decode(::Type{$A}, x::$U)
            if x > $(reinterpret(U, ub))
                throw(DecodeError{$A}(x))
            end
            return reinterpret($T, x)
        end
        decode(::Type{$A}, x::Unsigned) = decode($A, $U(x))
    end
end
