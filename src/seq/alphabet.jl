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

typealias NucleotideAlphabet Union{DNAAlphabet,RNAAlphabet}

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

alphabet(::Type{DNAAlphabet{2}}) = ACGT
alphabet(::Type{RNAAlphabet{2}}) = ACGU
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

EncodeError{A,T}(::Type{A}, val::T) = EncodeError{A,T}(val)

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

DecodeError{A,T}(::Type{A}, val::T) = DecodeError{A,T}(val)

function Base.showerror{A}(io::IO, err::DecodeError{A})
    print(io, "cannot decode ", err.val, " in ", A)
end


# DNA and RNA alphabets
# ---------------------

for A in (DNAAlphabet, RNAAlphabet)
    T = eltype(A)
    @eval begin
        # 2-bit encoding
        @inline function encode(::Type{$(A){2}}, nt::$(T))
            if count_ones(nt) != 1 || !isvalid(nt)
                throw(EncodeError($(A){2}, nt))
            end
            return convert(UInt8, trailing_zeros(nt))
        end
        @inline function decode(::Type{$(A){2}}, x::UInt8)
            if x > 0x03
                throw(DecodeError($(A){2}, x))
            end
            return reinterpret($(T), 0x01 << x)
        end
        @inline decode(::Type{$(A){2}}, x::Unsigned) = decode($(A){2}, UInt8(x))

        # 4-bit encoding
        @inline function encode(::Type{$(A){4}}, nt::$(T))
            if !isvalid(nt)
                throw(EncodeError($(A){4}, nt))
            end
            return reinterpret(UInt8, nt)
        end
        @inline function decode(::Type{$(A){4}}, x::UInt8)
            if !isvalid($(T), x)
                throw(DecodeError($(A){4}, x))
            end
            return reinterpret($(T), x)
        end
        @inline decode(::Type{$(A){4}}, x::Unsigned) = decode($(A){4}, UInt8(x))
    end
end


# AminoAcidAlphabet
# -----------------

@inline function encode(::Type{AminoAcidAlphabet}, aa::AminoAcid)
    if aa > AA_Gap
        throw(EncodeError(AminoAcidAlphabet, aa))
    end
    return reinterpret(UInt8, aa)
end

@inline function decode(::Type{AminoAcidAlphabet}, x::UInt8)
    if x > 0x1b
        throw(DecodeError(AminoAcidAlphabet, x))
    end
    return reinterpret(AminoAcid, x)
end

@inline function decode(::Type{AminoAcidAlphabet}, x::Unsigned)
    return decode(AminoAcidAlphabet, UInt8(x))
end


# CharAlphabet
# ------------

@inline function encode(::Type{CharAlphabet}, char::Char)
    if char > '\U10ffff'
        throw(EncodeError(CharAlphabet, char))
    end
    return reinterpret(UInt32, char)
end

@inline function decode(::Type{CharAlphabet}, x::UInt32)
    if x > 0x10ffff
        throw(DecodeError(CharAlphabet, x))
    end
    return reinterpret(Char, x)
end

@inline function decode(::Type{CharAlphabet}, x::Unsigned)
    return decode(CharAlphabet, UInt32(x))
end
