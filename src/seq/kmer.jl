# K-mer
# =====

# A Kmer is a sequence <= 32nt, without any 'N's, packed in a single 64 bit value.
#
# While NucleotideSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.
#
# Nucleotides are filled from MSBs to LSBs and right-aligned.
# For example, the memory layout of "TACG" is:
#   64-bit: 0b 00 00 … 00 11 00 01 10
#    4-mer:                T  A  C  G

bitstype 64 Kmer{T<:Nucleotide, K}

typealias DNAKmer{K} Kmer{DNANucleotide, K}
typealias RNAKmer{K} Kmer{RNANucleotide, K}
typealias Codon RNAKmer{3}


# Conversion
# ----------

function convert{K}(::Type{DNAKmer{K}}, x::UInt64)
    mask = ~UInt64(0) >> (64 - 2K)
    return box(DNAKmer{K}, unbox(UInt64, x & mask))
end
function convert{K}(::Type{RNAKmer{K}}, x::UInt64)
    mask = ~UInt64(0) >> (64 - 2K)
    return box(RNAKmer{K}, unbox(UInt64, x & mask))
end
convert{K}(::Type{UInt64}, x::DNAKmer{K}) = box(UInt64, unbox(DNAKmer{K}, x))
convert{K}(::Type{UInt64}, x::RNAKmer{K}) = box(UInt64, unbox(RNAKmer{K}, x))
convert{K}(::Type{DNAKmer{K}}, x::RNAKmer{K}) = reinterpret(DNAKmer{K}, x)
convert{K}(::Type{RNAKmer{K}}, x::DNAKmer{K}) = reinterpret(RNAKmer{K}, x)


"Convert a NucleotideSequence to a Kmer"
function convert{T, K}(::Type{Kmer{T, K}},
                       seq::Union{AbstractString,NucleotideSequence{T}})
    return make_kmer(Kmer{T,K}, seq)
end

# create a kmer from a sequence whose elements are convertible to a nucleotide
function make_kmer{T,K}(::Type{Kmer{T,K}},
                        seq::Union{AbstractString,NucleotideSequence,NTuple{K,T}})
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a NucleotideSequence of length $(length(seq))")

    x = UInt64(0)
    for c in seq
        nt = convert(T, c)
        if nt == nnucleotide(T)
            error("A K-mer may not contain an N in its sequence")
        end
        x = x << 2 | UInt64(nt)
    end

    return Kmer{T,K}(x)
end

make_kmer{K,T}(seq::NTuple{K,T}) = make_kmer(Kmer{T,K}, seq)

convert{T}(::Type{Kmer{T}}, seq::AbstractString) = convert(Kmer{T, length(seq)}, seq)
convert{T}(::Type{Kmer},    seq::NucleotideSequence{T}) = convert(Kmer{T, length(seq)}, seq)
convert{T}(::Type{Kmer{T}}, seq::NucleotideSequence{T}) = convert(Kmer{T, length(seq)}, seq)

convert{T, K}(::Type{AbstractString}, seq::Kmer{T, K}) = convert(AbstractString, [convert(Char, x) for x in seq])

function convert{T, K}(::Type{NucleotideSequence{T}}, x::Kmer{T, K})
    return convert(NucleotideSequence{T}, [T(nt) for nt in x])
end
convert{T, K}(::Type{NucleotideSequence}, x::Kmer{T, K}) = convert(NucleotideSequence{T}, x)


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

function =={T, K}(a::NucleotideSequence{T}, b::Kmer{T, K})
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

function =={T, K}(a::Kmer{T, K}, b::NucleotideSequence{T})
    return b == a
end

Base.hash(x::Kmer, h::UInt) = hash(UInt64(x), h)

function getindex{T, K}(x::Kmer{T, K}, i::Integer)
    if i < 1 || i > K
        throw(BoundsError())
    end
    return unsafe_getindex(x, i)
end

@inline function unsafe_getindex{T,K}(x::Kmer{T,K}, i::Integer)
    return convert(T, (UInt64(x) >> (2K - 2i)) & 0b11)
end

function show{T,K}(io::IO, x::Kmer{T,K})
    println(io, (T === DNANucleotide ? "DNA " : "RNA "), K, "-mer:")
    for i in 1:K
        write(io, convert(Char, x[i]))
    end
end

isless{T, K}(x::Kmer{T, K}, y::Kmer{T, K}) = isless(UInt64(x), UInt64(y))

length{T, K}(x::Kmer{T, K}) = K
endof(x::Kmer) = length(x)

# Iterating over nucleotides
start(x::Kmer) = 1
next(x::Kmer, i::Int) = unsafe_getindex(x, i), i + 1
done(x::Kmer, i::Int) = i > length(x)


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
reverse{T, K}(x::Kmer{T, K}) = Kmer{T,K}(nucrev(UInt64(x)) >> (64 - 2K))

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
mismatches{T, K}(a::Kmer{T, K}, b::Kmer{T, K}) = nucmismatches(UInt64(a), UInt64(b))

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

start(it::KmerNeighborIterator) = UInt64(0)
done(it::KmerNeighborIterator, i) = i == 4
length(::KmerNeighborIterator) = 4
function next{T, K}(it::KmerNeighborIterator{T, K}, i)
    return Kmer{T,K}((UInt64(it.x) << 2) | i), i + 1
end


"""
`NucleotideCounts(seq::Kmer)`

Constructs a NucleotideCounts object from a Kmer `seq`.
"""
function NucleotideCounts{T,K}(seq::Kmer{T, K})
    x         = convert(UInt64, seq)
    counts    = NucleotideCounts{T}()
    counts.a += count_a(x) - 32 + K # Take leading zeros into account
    counts.c += count_c(x)
    counts.g += count_g(x)
    counts.t += count_t(x)
    return counts
end
