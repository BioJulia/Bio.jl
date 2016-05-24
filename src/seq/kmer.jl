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
typealias Codon RNAKmer{3}


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
function make_kmer{T,K}(::Type{Kmer{T,K}},
                        seq::Union{AbstractString,BioSequence,NTuple{K,T}})
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a BioSequence of length $(length(seq))")

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


# Constructors
# ------------

"Construct a DNAKmer to an AbstractString"
dnakmer(seq::AbstractString) = convert(DNAKmer, seq)

"Construct a RNAKmer to an AbstractString"
rnakmer(seq::AbstractString) = convert(RNAKmer, seq)

"Construct a Kmer from a sequence of Nucleotides"
function kmer{T <: Nucleotide}(nts::T...)
    return make_kmer(nts)
end

"Construct a Kmer from a DNASequence"
function kmer(seq::DNASequence)
    return DNAKmer{length(seq)}(seq)
end

"Construct a Kmer from a RNASequence"
function kmer(seq::RNASequence)
    return RNAKmer{length(seq)}(seq)
end

# call kmer with @inline macro would reduce the performance significantly?
# Would the compiler inline even without @inline?
"Construct a DNAKmer from a DNASequence"
function dnakmer(seq::DNASequence)
    return DNAKmer{length(seq)}(seq)
end

"Construct a RNAKmer from a RNASequence"
function rnakmer(seq::RNASequence)
    return RNAKmer{length(seq)}(seq)
end


# Basic Functions
# ---------------

alphabet{k}(::Type{DNAKmer{k}}) = DNA_A:DNA_T
alphabet{k}(::Type{RNAKmer{k}}) = RNA_A:RNA_U

for (A, N) in ((DNAAlphabet, DNANucleotide), (RNAAlphabet, RNANucleotide))
    @compat @eval begin
        function Base.:(==){A<:$A,K}(a::BioSequence{A}, b::Kmer{$N,K})
            if length(a) != K
                return false
            end
            for (u, v) in zip(a, b)
                if u != v
                    return false
                end
            end
            return true
        end
        Base.:(==){A<:$A,K}(a::Kmer{$N,K}, b::BioSequence{A}) = b == a
    end
end

Base.hash(x::Kmer, h::UInt) = hash(UInt64(x), h)

function Base.checkbounds(x::Kmer, i::Integer)
    if 1 ≤ i ≤ endof(x)
        return true
    end
    throw(BoundsError(x, i))
end

function Base.getindex{T, K}(x::Kmer{T, K}, i::Integer)
    checkbounds(x, i)
    return unsafe_getindex(x, i)
end

@inline function unsafe_getindex{T,K}(x::Kmer{T,K}, i::Integer)
    return convert(T, (UInt64(x) >> (2K - 2i)) & 0b11)
end

Base.summary{k}(x::DNAKmer{k}) = string("DNA ", k, "-mer")
Base.summary{k}(x::RNAKmer{k}) = string("RNA ", k, "-mer")

function Base.show(io::IO, x::Kmer)
    println(io, summary(x), ':')
    print(io, x)
end

Base.print(io::IO, x::Kmer) = showcompact(io, x)

function Base.showcompact{T,k}(io::IO, x::Kmer{T,k})
    for i in 1:k
        write(io, convert(Char, x[i]))
    end
end

@compat begin
    Base.:-{T,K}(x::Kmer{T,K}, y::Integer) =
        Kmer{T,K}(UInt64(x) - reinterpret(UInt64, Int64(y)))
    Base.:+{T,K}(x::Kmer{T,K}, y::Integer) =
        Kmer{T,K}(UInt64(x) + reinterpret(UInt64, Int64(y)))
    Base.:+{T,K}(x::Integer, y::Kmer{T,K}) = y + x
end
Base.isless{T,K}(x::Kmer{T,K}, y::Kmer{T,K}) = isless(UInt64(x), UInt64(y))

Base.length{T, K}(x::Kmer{T, K}) = K
Base.endof(x::Kmer) = length(x)

# Iterating over nucleotides
Base.start(x::Kmer) = 1
Base.next(x::Kmer, i::Int) = unsafe_getindex(x, i), i + 1
Base.done(x::Kmer, i::Int) = i > length(x)
Base.eltype{T,k}(::Type{Kmer{T,k}}) = T

function Base.findnext{T}(kmer::Kmer{T}, val, start::Integer)
    checkbounds(kmer, start)
    v = convert(T, val)
    for i in Int(start):endof(kmer)
        x = unsafe_getindex(kmer, i)
        if x == v
            return i
        end
    end
    return 0
end

function Base.findprev{T}(kmer::Kmer{T}, val, start::Integer)
    checkbounds(kmer, start)
    v = convert(T, val)
    for i in Int(start):-1:1
        x = unsafe_getindex(kmer, i)
        if x == v
            return i
        end
    end
    return 0
end

Base.findfirst(kmer::Kmer, val) = findnext(kmer, val, 1)
Base.findlast(kmer::Kmer, val)  = findprev(kmer, val, endof(kmer))


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

Base.start(it::KmerNeighborIterator) = UInt64(0)
Base.done(it::KmerNeighborIterator, i) = i == 4
Base.length(::KmerNeighborIterator) = 4
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
