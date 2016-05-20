# K-mer
# =====
#
# Compact k-mer sequence type.
#
# A Kmer is a sequence ≤ 32nt, without any 'N's, packed in a single 64 bit
# value.  While BioSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Representation
# --------------
#
# Four kinds of nucleotides are encoded as follows:
#
#   nucleotide | binary
#   ---------- | ------
#       A      |   00
#       C      |   01
#       G      |   10
#     T / U    |   11
#
# Nucleotides are filled from MSBs to LSBs and right-aligned so that all k-mers
# are lexicographically ordered. For example, the memory layout of "TACG" is:
#   64-bit: 0b 00 00 … 00 11 00 01 10
#    4-mer:                T  A  C  G

bitstype 64 Kmer{T<:Nucleotide, K} <: Sequence

typealias DNAKmer{K} Kmer{DNANucleotide, K}
typealias RNAKmer{K} Kmer{RNANucleotide, K}
typealias DNACodon DNAKmer{3}
typealias RNACodon RNAKmer{3}

function Kmer{T<:Nucleotide}(nts::T...)
    return make_kmer(nts)
end

function DNACodon(x::DNANucleotide, y::DNANucleotide, z::DNANucleotide)
    return make_kmer((x, y, z))
end

function RNACodon(x::RNANucleotide, y::RNANucleotide, z::RNANucleotide)
    return make_kmer((x, y, z))
end


# Conversion
# ----------

function Base.convert{K}(::Type{DNAKmer{K}}, x::UInt64)
    mask = ~UInt64(0) >> (64 - 2K)
    return reinterpret(DNAKmer{K}, x & mask)
end
function Base.convert{K}(::Type{RNAKmer{K}}, x::UInt64)
    mask = ~UInt64(0) >> (64 - 2K)
    return reinterpret(RNAKmer{K}, x & mask)
end
Base.convert{K}(::Type{UInt64}, x::DNAKmer{K}) = reinterpret(UInt64, x)
Base.convert{K}(::Type{UInt64}, x::RNAKmer{K}) = reinterpret(UInt64, x)
Base.convert{K}(::Type{DNAKmer{K}}, x::RNAKmer{K}) = reinterpret(DNAKmer{K}, x)
Base.convert{K}(::Type{RNAKmer{K}}, x::DNAKmer{K}) = reinterpret(RNAKmer{K}, x)

function Base.convert{T,K}(::Type{Kmer{T,K}}, seq::AbstractString)
    return make_kmer(Kmer{T,K}, seq)
end

function Base.convert{T,K,A<:DNAAlphabet}(::Type{Kmer{T,K}}, seq::BioSequence{A})
    return make_kmer(Kmer{DNANucleotide,K}, seq)
end

function Base.convert{T,K,A<:RNAAlphabet}(::Type{Kmer{T,K}}, seq::BioSequence{A})
    return make_kmer(Kmer{RNANucleotide,K}, seq)
end

Base.convert{T}(::Type{Kmer{T}}, seq::AbstractString) = convert(Kmer{T,length(seq)}, seq)
Base.convert{A<:DNAAlphabet}(::Type{Kmer}, seq::BioSequence{A}) = convert(Kmer{DNANucleotide,length(seq)}, seq)
Base.convert{A<:RNAAlphabet}(::Type{Kmer}, seq::BioSequence{A}) = convert(Kmer{RNANucleotide,length(seq)}, seq)
Base.convert{A<:DNAAlphabet}(::Type{DNAKmer}, seq::BioSequence{A}) = convert(DNAKmer{length(seq)}, seq)
Base.convert{A<:RNAAlphabet}(::Type{RNAKmer}, seq::BioSequence{A}) = convert(RNAKmer{length(seq)}, seq)

# create a kmer from a sequence whose elements are convertible to a nucleotide
function make_kmer{T,K}(::Type{Kmer{T,K}}, seq)
    seqlen = length(seq)
    if seqlen > 32
        throw(ArgumentError("cannot create a k-mer loger than 32nt"))
    elseif seqlen != K
        throw(ArgumentError("cannot create a $(K)-mer from a sequence of length $(seqlen)"))
    end

    x = UInt64(0)
    for c in seq
        nt = convert(T, c)
        if isambiguous(nt)
            error("A K-mer cannot contain an ambiguous nucleotide in its sequence")
        end
        x = x << 2 | UInt64(nt)
    end
    return Kmer{T,K}(x)
end

make_kmer{K,T}(seq::NTuple{K,T}) = make_kmer(Kmer{T,K}, seq)

Base.convert{K}(::Type{BioSequence}, x::DNAKmer{K}) = DNASequence(x)
Base.convert{K}(::Type{BioSequence}, x::RNAKmer{K}) = RNASequence(x)
Base.convert{A<:DNAAlphabet,K}(::Type{BioSequence{A}}, x::DNAKmer{K}) = BioSequence{A}([nt for nt in x])
Base.convert{A<:RNAAlphabet,K}(::Type{BioSequence{A}}, x::RNAKmer{K}) = BioSequence{A}([nt for nt in x])
Base.convert{S<:AbstractString}(::Type{S}, seq::Kmer) = convert(S, [Char(x) for x in seq])


# Basic Functions
# ---------------

alphabet{k}(::Type{DNAKmer{k}}) = DNA_A:DNA_T
alphabet{k}(::Type{RNAKmer{k}}) = RNA_A:RNA_U

Base.hash(x::Kmer, h::UInt) = hash(UInt64(x), h)

Base.length{T,K}(x::Kmer{T, K}) = K
Base.eltype{T,k}(::Type{Kmer{T,k}}) = T

@inline function inbounds_getindex{T,K}(x::Kmer{T,K}, i::Integer)
    return convert(T, (UInt64(x) >> (2K - 2i)) & 0b11)
end

Base.summary{k}(x::DNAKmer{k}) = string("DNA ", k, "-mer")
Base.summary{k}(x::RNAKmer{k}) = string("RNA ", k, "-mer")

@compat begin
    Base.:-{T,K}(x::Kmer{T,K}, y::Integer) =
        Kmer{T,K}(UInt64(x) - reinterpret(UInt64, Int64(y)))
    Base.:+{T,K}(x::Kmer{T,K}, y::Integer) =
        Kmer{T,K}(UInt64(x) + reinterpret(UInt64, Int64(y)))
    Base.:+{T,K}(x::Integer, y::Kmer{T,K}) = y + x
    Base.:(==){T,k}(x::Kmer{T,k}, y::Kmer{T,k}) = UInt64(x) == UInt64(y)
end
Base.isless{T,K}(x::Kmer{T,K}, y::Kmer{T,K}) = isless(UInt64(x), UInt64(y))


# Other functions
# ---------------

"""
`complement(kmer::Kmer)`

The Kmer complement of `kmer`
"""
complement{T, K}(x::Kmer{T, K}) = Kmer{T,K}(~UInt64(x))

"""
`reverse(kmer::Kmer)`

Reversed copy of `kmer`
"""
Base.reverse{T, K}(x::Kmer{T, K}) = Kmer{T,K}(nucrev2(UInt64(x)) >> (64 - 2K))

"""
`reverse_complement(kmer::Kmer)`

Reversed complement of `kmer`
"""
reverse_complement{T, K}(x::Kmer{T, K}) = complement(reverse(x))

"""
`mismatches(a::Kmer, b::Kmer)`

Return the number of mismatches between `a` and `b`.

### Arguments
* `a`: first sequence to compare
* `b`: second sequence to compare

### Returns
The number of mismatches
"""
mismatches{T, K}(a::Kmer{T, K}, b::Kmer{T, K}) = nuc2mismatches(UInt64(a), UInt64(b))

"""
Canonical k-mer of `x`

A canonical k-mer is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting k-mers in data that is not strand specific,
and thus observing k-mer is equivalent to observing its reverse complement.
"""
function canonical{T, K}(x::Kmer{T, K})
    y = reverse_complement(x)
    return x < y ? x : y
end


# De Bruijn Graph
# ---------------

# neighbors on a de Bruijn graph
immutable KmerNeighborIterator{T, K}
    x::Kmer{T, K}
end

"""
Iterate through k-mers neighboring `x` on a de Bruijn graph.
"""
function neighbors(x::Kmer)
    return KmerNeighborIterator(x)
end

Base.length(::KmerNeighborIterator) = 4
Base.eltype{T,k}(::Type{KmerNeighborIterator{T,k}}) = Kmer{T,k}
Base.start(it::KmerNeighborIterator) = UInt64(0)
Base.done(it::KmerNeighborIterator, i) = i == 4
function Base.next{T, K}(it::KmerNeighborIterator{T, K}, i)
    return Kmer{T,K}((UInt64(it.x) << 2) | i), i + 1
end

# Count A, C, T/U, G respectively in a kmer stored in a UInt64
function count_a(x::UInt64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & 0x5555555555555555)
end
count_c(x::UInt64) = count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
count_g(x::UInt64) = count_ones(((x >>> 1) & (~x)) & 0x5555555555555555)
count_t(x::UInt64) = count_ones((x    & (x >>> 1)) & 0x5555555555555555)
