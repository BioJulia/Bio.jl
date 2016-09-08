# Nucleotides
# ===========
#
# DNA and RNA nucleotide types.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Nucleotide Encoding
# -------------------
#
# Unambiguous nucleotides are represented in one-hot encoding as follows:
#
#   | Nucleotide | Bits |
#   | ---------- | ---- |
#   |     A      | 0001 |
#   |     C      | 0010 |
#   |     G      | 0100 |
#   |    T/U     | 1000 |
#
# Ambiguous nucleotides are bitwise OR of these four nucleotides. For example, R
# , A or G, is represented as 0101 (= A: 0001 | G: 0100). The gap symbol is
# always 0000.  The meaningful four bits are stored in the least significant
# bits of a byte.

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


# Bit Operations
# --------------

function Base.:~{N<:Nucleotide}(x::N)
    return reinterpret(N, ~reinterpret(UInt8, x) & 0b1111)
end

function Base.:|{N<:Nucleotide}(x::N, y::N)
    return reinterpret(N, reinterpret(UInt8, x) | reinterpret(UInt8, y))
end

function Base.:&{N<:Nucleotide}(x::N, y::N)
    return reinterpret(N, reinterpret(UInt8, x) & reinterpret(UInt8, y))
end

function Base.:-{N<:Nucleotide}(x::N, y::N)
    return Int(x) - Int(y)
end

function Base.:-{N<:Nucleotide}(x::N, y::Integer)
    return x + (-y)
end

function Base.:+{N<:Nucleotide}(x::N, y::Integer)
    return reinterpret(N, (UInt8(x) + y % UInt8) & 0b1111)
end

function Base.isless{N<:Nucleotide}(x::N, y::N)
    return isless(reinterpret(UInt8, x), reinterpret(UInt8, y))
end

@inline function Base.count_ones(nt::Nucleotide)
    return count_ones(reinterpret(UInt8, nt))
end

function Base.trailing_zeros(nt::Nucleotide)
    return trailing_zeros(reinterpret(UInt8, nt))
end

function gap{N<:Nucleotide}(::Type{N})
    return reinterpret(N, 0b0000)
end

"""
    isGC(nt::Nucleotide)

Test if `nt` is surely either guanine or cytosine.
"""
function isGC(nt::Nucleotide)
    bits = reinterpret(UInt8, nt)
    return bits != 0 && (bits & 0b1001) == 0
end

"""
    ispurine(nt::Nucleotide)

Test if nucleotide is surely a purine.
"""
@inline function ispurine(nt::Nucleotide)
    bits = reinterpret(UInt8, nt)
    return bits != 0 && (bits & 0b1010) == 0
end

"""
    ispyrimidine(nt::Nucleotide)

Test if nucleotide is surely a pyrimidine.
"""
@inline function ispyrimidine(nt::Nucleotide)
    bits = reinterpret(UInt8, nt)
    return bits != 0 && (bits & 0b0101) == 0
end

"""
    isambiguous(nt::Nucleotide)

Test if `nt` is ambiguous nucleotide.
"""
@inline function isambiguous(nt::Nucleotide)
    return count_ones(nt) > 1
end

"""
    complement(nt::Nucleotide)

Return the complementary nucleotide of `nt`.
"""
function complement(nt::Nucleotide)
    bits = compatbits(nt)
    return reinterpret(
        typeof(nt),
        (bits & 0x01) << 3 | (bits & 0x08) >> 3 |
        (bits & 0x02) << 1 | (bits & 0x04) >> 1)
end

function Base.isvalid{T<:Nucleotide}(::Type{T}, x::Integer)
    return 0 ≤ x < 16
end

function Base.isvalid(nt::Nucleotide)
    return reinterpret(UInt8, nt) ≤ 0b1111
end

"""
    iscompatible(x, y)

Return `true` if and only if `x` and `y` are compatible with each other (i.e.
`x` and `y` can be the same symbol).

`x` and `y` must be the same type (`DNANucleotide`, `RNANucleotide` or `AminoAcid`).

# Examples

```julia
julia> iscompatible(DNA_A, DNA_A)
true

julia> iscompatible(DNA_C, DNA_N)  # DNA_N can be DNA_C
true

julia> iscompatible(DNA_C, DNA_R)  # DNA_R (A or G) cannot be DNA_C
false

julia> iscompatible(AA_A, AA_X)    # AA_X can be AA_A
true

```
"""
@inline function iscompatible{T<:Nucleotide}(x::T, y::T)
    return compatbits(x) & compatbits(y) != 0
end

# Return the compatibility bits of `nt`.
@inline function compatbits(nt::Nucleotide)
    return reinterpret(UInt8, nt)
end


# Nucleotide encoding definition
# ------------------------------

# DNA Nucleotides

# lookup table for characters
const char_to_dna = [0x80 for _ in 0x00:0xff]
const dna_to_char = Vector{Char}(16)

# derived from "The DDBJ/ENA/GenBank Feature Table Definition"
# §7.4.1 Nucleotide base code (IUPAC)
# http://www.insdc.org/documents/feature_table.html#7.4.1
for (char, doc, bits) in [
        ('-', "DNA Gap",                                   0b0000),
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
        ('N', "DNA Adenine, Cytosine, Guanine or Thymine", 0b1111)]
    var = Symbol("DNA_", char != '-' ? char : "Gap")
    @eval begin
        @doc $(doc) const $(var) = reinterpret(DNANucleotide, $(bits))
        char_to_dna[$(Int(char)+1)] = char_to_dna[$(Int(lowercase(char))+1)] = $(bits)
        dna_to_char[$(Int(bits)+1)] = $(char)
    end
end

@eval alphabet(::Type{DNANucleotide}) = $(tuple([reinterpret(DNANucleotide, x)
                                                 for x in 0b0000:0b1111]...))

const ACGT = (DNA_A, DNA_C, DNA_G, DNA_T)
const ACGTN = (DNA_A, DNA_C, DNA_G, DNA_T, DNA_N)

# RNA Nucleotides

# lookup table for characters
const char_to_rna = [0x80 for _ in 0x00:0xff]
const rna_to_char = Vector{Char}(16)

for (char, doc, dna) in [
        ('-', "RNA Gap",                                  DNA_Gap),
        ('A', "RNA Adenine",                              DNA_A  ),
        ('C', "RNA Cytosine",                             DNA_C  ),
        ('G', "RNA Guanine",                              DNA_G  ),
        ('U', "RNA Uracil",                               DNA_T  ),
        ('M', "RNA Adenine or Cytosine",                  DNA_M  ),
        ('R', "RNA Adenine or Guanine",                   DNA_R  ),
        ('W', "RNA Adenine or Uracil",                    DNA_W  ),
        ('S', "RNA Cytosine or Guanine",                  DNA_S  ),
        ('Y', "RNA Cytosine or Uracil",                   DNA_Y  ),
        ('K', "RNA Guanine or Uracil",                    DNA_K  ),
        ('V', "RNA Adenine, Cytosine or Guanine",         DNA_V  ),
        ('H', "RNA Adenine, Cytosine or Uracil",          DNA_H  ),
        ('D', "RNA Adenine, Guanine or Uracil",           DNA_D  ),
        ('B', "RNA Cytosine, Guanine or Uracil",          DNA_B  ),
        ('N', "RNA Adenine, Cytosine, Guanine or Uracil", DNA_N  )]
    var = Symbol("RNA_", char != '-' ? char : "Gap")
    @eval begin
        @doc $(doc) const $(var) = reinterpret(RNANucleotide, $(dna))
        char_to_rna[$(Int(char)+1)] = char_to_rna[$(Int(lowercase(char)+1))] = reinterpret(UInt8, $(dna))
        rna_to_char[$(Int(dna)+1)] = $(char)
    end
end

@eval alphabet(::Type{RNANucleotide}) = $(tuple([reinterpret(RNANucleotide, x)
                                                 for x in 0b0000:0b1111]...))

const ACGU = (RNA_A, RNA_C, RNA_G, RNA_U)
const ACGUN = (RNA_A, RNA_C, RNA_G, RNA_U, RNA_N)


# Print functions
# ---------------

function Base.convert(::Type{DNANucleotide}, c::Char)
    if c > '\xff'
        throw(InexactError())
    end
    @inbounds dna = char_to_dna[Int(c) + 1]
    if !isvalid(DNANucleotide, dna)
        throw(InexactError())
    end
    return reinterpret(DNANucleotide, dna)
end

function Base.convert(::Type{RNANucleotide}, c::Char)
    if c > '\xff'
        throw(InexactError())
    end
    @inbounds rna = char_to_rna[Int(c) + 1]
    if !isvalid(RNANucleotide, rna)
        throw(InexactError())
    end
    return reinterpret(RNANucleotide, rna)
end

function Base.convert(::Type{Char}, nt::DNANucleotide)
    return dna_to_char[convert(UInt8, nt) + 1]
end

function Base.convert(::Type{Char}, nt::RNANucleotide)
    return rna_to_char[convert(UInt8, nt) + 1]
end

function Base.show(io::IO, nt::DNANucleotide)
    if isvalid(nt)
        if nt == DNA_Gap
            write(io, "DNA_Gap")
        else
            write(io, "DNA_", Char(nt))
        end
    else
        write(io, "Invalid DNA Nucleotide")
    end
    return
end

function Base.show(io::IO, nt::RNANucleotide)
    if isvalid(nt)
        if nt == RNA_Gap
            write(io, "RNA_Gap")
        else
            write(io, "RNA_", Char(nt))
        end
    else
        write(io, "Invalid RNA Nucleotide")
    end
    return
end

function Base.print(io::IO, nt::Nucleotide)
    if !isvalid(nt)
        throw(ArgumentError("nucleotide is invalid"))
    end
    write(io, Char(nt))
    return
end
