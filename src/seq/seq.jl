
module Seq

using Base.Intrinsics
using IntervalTrees

import Base: convert, getindex, show, length, start, next, done, copy, reverse

export Nucleotide, DNANucleotide, RNANucleotide,
       NucleotideSequence, DNASequence, RNASequence, @dna_str, @rna_str,
       complement, reverse_complement


# Nucleotides
# -----------


abstract Nucleotide
bitstype 8 DNANucleotide <: Nucleotide
bitstype 8 RNANucleotide <: Nucleotide


function convert(::Type{DNANucleotide}, nt::Uint8)
    return box(DNANucleotide, unbox(Uint8, nt))
end


function convert(::Type{RNANucleotide}, nt::Uint8)
    return box(RNANucleotide, unbox(Uint8, nt))
end


function convert(::Type{Uint8}, nt::DNANucleotide)
    return box(Uint8, unbox(DNANucleotide, nt))
end


function convert(::Type{Uint8}, nt::RNANucleotide)
    return box(Uint8, unbox(RNANucleotide, nt))
end


function convert{T <: Unsigned, S <: Nucleotide}(::Type{T}, nt::S)
    return box(T, Base.zext_int(T, unbox(S, nt)))
end



# TODO: Should we use ape-style bit masks for these?
const DNA_A = convert(DNANucleotide, 0b000)
const DNA_C = convert(DNANucleotide, 0b001)
const DNA_G = convert(DNANucleotide, 0b010)
const DNA_T = convert(DNANucleotide, 0b011)
const DNA_N = convert(DNANucleotide, 0b100)
const DNA_INVALID = convert(DNANucleotide, 0b1000)

const RNA_A = convert(RNANucleotide, 0b000)
const RNA_C = convert(RNANucleotide, 0b001)
const RNA_G = convert(RNANucleotide, 0b010)
const RNA_U = convert(RNANucleotide, 0b011)
const RNA_N = convert(RNANucleotide, 0b100)
const RNA_INVALID = convert(DNANucleotide, 0b1000)


function nnucleotide(::Type{DNANucleotide})
    return DNA_N
end


function nnucleotide(::Type{RNANucleotide})
    return RNA_N
end


const dna_to_char = ['A', 'C', 'G', 'T', 'N']

function convert(::Type{Char}, nt::DNANucleotide)
    return dna_to_char[convert(Uint8, nt) + 1]
end

# lookup table for characters in 'A':'n'
const char_to_dna = [
     DNA_A,       DNA_INVALID, DNA_C,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_N,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_T,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_A,       DNA_INVALID, DNA_C,       DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_N ]

function convert(::Type{DNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 'n' ? char_to_dna[c - 'A' + 1] : DNA_INVALID
    if nt == DNA_INVALID
        error("$(c) is not a valid DNA nucleotide")
    end

    return nt
end


const rna_to_char = ['A', 'C', 'G', 'U', 'N']

function convert(::Type{Char}, nt::RNANucleotide)
    return rna_to_char[convert(Uint8, nt) + 1]
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

function convert(::Type{RNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 'u' ? char_to_rna[c - 'A' + 1] : RNA_INVALID
    if nt == RNA_INVALID
        error("$(c) is not a valid DNA nucleotide")
    end

    return nt
end


function show(io::IO, nt::DNANucleotide)
    write(io, convert(Char, nt))
end


function show(io::IO, nt::RNANucleotide)
    write(io, convert(Char, nt))
end


# Nucleotide Sequences
# --------------------

# How many Uint64s are needed to represent a sequence of length n
function seq_data_len(n::Integer)
    d, r = divrem(n, 32)
    return d + (r > 0 ? 1 : 0)
end

# A representation of DNA and RNA sequences.
#
# NucleotideSequence uses 2-bit encoding to efficiently store and operate on
# sequences. Despite the encoding, N's are still allowed in the sequence and
# represented by an N “mask” stored in an IntervalTree. This has the following
# tradeoff: if Ns are relatively rare or occur in large blocks (as in most
# genome sequences), this is an extremely space efficient representation.
# Performance degrades in long sequences with many interspersed Ns, but we hope
# to an acceptable degree.
#
# NucleotideSequence is effectively immutable. Though it's possible to do so, we
# do not export functions that modify sequences in place. Immutability allows
# for efficient subsequence operations. Every NucleotideSequence has a `part`
# field that specifies a range within the underlying data. When subsequences are
# made, rather than coping the data, a new `part` is specified.
#
immutable NucleotideSequence{T <: Nucleotide}
    # 2-bit encoded sequence
    data::Vector{Uint64}

    # 'N' mask
    ns::BitVector

    # interval within data defining the (sub)sequence
    part::UnitRange{Int}

    # Construct from raw components
    function NucleotideSequence(data::Vector{Uint64}, ns::BitVector,
                                part::UnitRange)
        return new(data, ns, part)
    end

    # Construct a subsequence of another nucleotide sequence
    function NucleotideSequence(other::NucleotideSequence, part::UnitRange)
        start = other.part.start + part.start - 1
        stop = start + length(part) - 1
        @assert start >= other.part.start
        @assert stop <= other.part.stop
        return new(other.data, other.ns, part)
    end

    # Construct an empty sequence
    function NucleotideSequence()
        return new(zeros(Uint64, 0), BitVector(0), 1:0)
    end

    # Construct a sequence from a string
    function NucleotideSequence(seq::String)
        len = seq_data_len(length(seq))
        data = zeros(Uint64, len)
        ns = BitArray(length(seq))
        fill!(ns, false)

        j = start(seq)
        idx = 1
        @inbounds for i in 1:length(data)
            shift = 0
            while shift < 64 && !done(seq, j)
                c, j = next(seq, j)
                nt = convert(T, c)
                if nt == nnucleotide(T)
                    ns[idx] = true
                else
                    data[i] |= convert(Uint64, nt) << shift
                end
                idx += 1
                shift += 2
            end
        end

        return new(data, ns, 1:length(seq))
    end
end


typealias DNASequence NucleotideSequence{DNANucleotide}
typealias RNASequence NucleotideSequence{RNANucleotide}

# Copy a sequence.
#
# Unlike constructing subsequences with seq[a:b], this function actually copies
# the underlying data. Since sequences are immutable, you should basically
# never have to do this. It's useful only for low-level algorithms like
# `reverse` which actually do make a copy and modify the copy in
# place.
function copy{T}(seq::NucleotideSequence{T})
    data = zeros(Uint64, seq_data_len(length(seq)))
    d0, r0 = divrem(seq.part.start - 1, 32)

    h = 64 - 2*r0
    k = 2*r0

    j = d0 + 1
    for i in 1:length(data)
        data[i] |= seq.data[j] >>> k

        j += 1
        if j > length(seq.data)
            break
        end

        data[i] |= seq.data[j] << h
    end

    return NucleotideSequence{T}(data, seq.ns[seq.part.start:seq.part.stop],
                                 1:length(seq))
end


# Iterating throug nucleotide sequences
function start(seq::NucleotideSequence)
    i = seq.part.start - 1
    return i
end


function next{T}(seq::NucleotideSequence{T}, i)
    nvalue, _ = next(seq.ns, i)
    if nvalue
        return (nnucleotide(T), i + 1)
    else
        return (getnuc(T, seq.data, i + 1), i + 1)
    end
end


function done(seq::NucleotideSequence, i)
    return i >= length(seq)
end


# String decorator syntax to enable building sequence literals like:
#     dna"ACGTACGT" and rna"ACGUACGU"
macro dna_str(seq, flags...)
    return DNASequence(seq)
end

macro rna_str(seq, flags...)
    return RNASequence(seq)
end


function length(seq::NucleotideSequence)
    return length(seq.part)
end


function show{T}(io::IO, seq::NucleotideSequence{T})
    if T == DNANucleotide
        write(io, "dna\"")
    elseif T == RNANucleotide
        write(io, "rna\"")
    end

    len = length(seq)
    for nt in seq
        write(io, convert(Char, nt))
    end

    write(io, "\"")
end


# get nucleotide at position i, ignoring the N mask
function getnuc(T::Type, data::Vector{Uint64}, i::Integer)
    d, r = divrem(i - 1, 32)
    return convert(T, convert(Uint8, (data[d + 1] >>> (2*r)) & 0b11))
end


function getindex{T}(seq::NucleotideSequence{T}, i::Integer)
    if i > length(seq) || i < 1
        error(BoundsError())
    end
    i += seq.part.start - 1
    if seq.ns[i]
        return nnucleotide(T)
    else
        return getnuc(T, seq.data, i)
    end
end


function getindex{T}(seq::NucleotideSequence{T}, r::UnitRange)
    return NucleotideSequence{T}(seq, r)
end


function convert(::Type{DNASequence}, seq::String)
    return DNASequence(seq)
end


function convert(::Type{RNASequence}, seq::String)
    return DNASequence(seq)
end


function convert(::Type{String}, seq::NucleotideSequence)
    return convert(String, [convert(Char, x) for x in seq])
end


# Transformations
# ---------------


# In-place complement for low level work.
function _complement!(seq::NucleotideSequence)
    for i in 1:length(seq.data)
        seq.data[i] = ~seq.data[i]
    end
end


function complement(seq::NucleotideSequence)
    seq = copy(seq)
    _complement!(seq)
    return seq
end


# Nucleotide reverse. Reverse a kmer stored in a Uint64.
function nucrev(x::Uint64)
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>>  2
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >>> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >>> 32
     return x
end


function reverse{T}(seq::NucleotideSequence{T})
    k = (2 * length(seq)) % 64
    h = 64 - k
    if k == 0
        k = 64
        h = 0
    end

    data = zeros(Uint64, length(seq.data))
    j = length(data)
    for i in 1:length(data)
        x = nucrev(seq.data[i])
        data[j] |= x >>> h
        if (j -= 1) == 0
            break
        end
        data[j] |= x << k;
    end

    return NucleotideSequence{T}(data, reverse(seq.ns), seq.part)
end


function reverse_complement(seq::NucleotideSequence)
    seq = reverse(seq)
    _complement!(seq)
    return seq
end


end
