
module Seq

using Base.Intrinsics

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof

export Nucleotide, DNANucleotide, RNANucleotide,
       NucleotideSequence, DNASequence, RNASequence, @dna_str, @rna_str,
       complement, reverse_complement, mismatches, ns, eachsubseq, eachkmer


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
type NucleotideSequence{T <: Nucleotide}
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

# If a sequence's starting position within its data != 1, this function
# will copy a subset of the data to align with sequences range and make start == 1.
#
# The user should never need to call this, as it has no outward effect on the
# sequence, but it makes functions like mismatch easier and faster if can assume
# a sequence is aligned with its data.
#
# If reorphan = true, copy the data regardless of the start position.
#
function orphan!{T}(seq::NucleotideSequence{T}, reorphan=false)
    if !reorphan && seq.part.start == 1
        return
    end

    data = zeros(Uint64, seq_data_len(length(seq)))
    d0, r0 = divrem(seq.part.start - 1, 32)

    h = 64 - 2*r0
    k = 2*r0

    j = d0 + 1
    @inbounds for i in 1:length(data)
        data[i] |= seq.data[j] >>> k

        j += 1
        if j > length(seq.data)
            break
        end

        data[i] |= seq.data[j] << h
    end

    seq.data = data
    seq.ns = seq.ns[seq.part.start:seq.part.stop]
    seq.part = 1:length(seq.part)
    return seq
end

# Copy a sequence.
#
# Unlike constructing subsequences with seq[a:b], this function actually copies
# the underlying data. Since sequences are immutable, you should basically
# never have to do this. It's useful only for low-level algorithms like
# `reverse` which actually do make a copy and modify the copy in
# place.
function copy{T}(seq::NucleotideSequence{T})
    return orphan!(NucleotideSequence{T}(seq.data, seq.ns, seq.part), true)
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
    return i >= seq.part.stop
end


# Iterate through positions of Ns efficiently
# TODO: This is code that should ultimately be put in the bitarray code Base.

immutable SequenceNIterator
    ns::BitVector
    part::UnitRange{Int}
end


function ns(seq::NucleotideSequence)
    return SequenceNIterator(seq.ns, seq.part)
end


function nextn(it::SequenceNIterator, i)
    d, r = divrem(i - 1, 64)
    while d < length(it.ns.chunks) && it.ns.chunks[d + 1] >>> r == 0 && d * 64 < it.part.stop
        d += 1
        r = 0
    end

    if d * 64 + r + 1 > it.part.stop
        return d * 64 + r + 1
    end

    if d + 1 <= length(it.ns.chunks)
        x = it.ns.chunks[d + 1] >>> r
        while x & 0x1 == 0
            x >>>= 1
            r += 1
        end
    end

    return d * 64 + r + 1
end


function start(it::SequenceNIterator)
    return nextn(it, it.part.start)
end


function next(it::SequenceNIterator, i)
    d, r = divrem(i - 1, 64)
    next_i = nextn(it, i + 1)
    return i + it.part.start - 1, next_i
end


function done(it::SequenceNIterator, i)
    return i > it.part.stop
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


function endof(seq::NucleotideSequence)
    return length(seq)
end


function show{T}(io::IO, seq::NucleotideSequence{T})

    # don't show more than this many characters to avoid
    const maxcount = 60
    if T == DNANucleotide
        write(io, "dna\"")
    elseif T == RNANucleotide
        write(io, "rna\"")
    end

    len = length(seq)
    if len > maxcount
        for nt in seq[1:div(maxcount, 2) - 1]
            write(io, convert(Char, nt))
        end
        write("…")
        for nt in seq[(end - (div(maxcount, 2) - 1)):end]
            write(io, convert(Char, nt))
        end
    else
        for nt in seq
            write(io, convert(Char, nt))
        end
    end

    write(io, "\"  # ", string(len), "nt ",
          T == DNANucleotide ? "DNA" : "RNA", " sequence")

end


# get nucleotide at position i, ignoring the N mask
function getnuc(T::Type, data::Vector{Uint64}, i::Integer)
    d, r = divrem(i - 1, 32)
    return convert(T, convert(Uint8, (data[d + 1] >>> (2*r)) & 0b11))
end


function zeronuc!(T::Type, data::Vector{Uint64}, i::Integer)
    d, r = divrem(i - 1, 32)
    data[d + 1] &= ~(convert(Uint64, nt) << (2*r))
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


# TODO: do we actually need this
# Convert a sequence to an integer if possible.
#
# This is slightly odd, but ends up being quite useful for indexing into arrays,
# or counting k-mers.
#
#function convert{T <: Unsigned}(::Type{T}, seq::NucleotideSequence)
    #if length(seq) > 4 * sizeof(T)
        #error("nucleotide sequence is too long to be converted to integer type $(T)")
    #end

    #if hasn(seq)
        #error("cannot convert a sequence with Ns to integer type $(T)")
    #end

    #x = zero(T)
    #for nt in seq
        #x = (x << 2) | convert(Uint8, nt)
    #end
    #return x
#end



# Transformations
# ---------------


# In-place complement for low level work.
function _complement!(seq::NucleotideSequence)
    @inbounds for i in 1:length(seq.data)
        seq.data[i] = ~seq.data[i]
    end
    # TODO: zero nucleotides with n
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
    orphan!(seq)

    k = (2 * length(seq)) % 64
    h = 64 - k
    if k == 0
        k = 64
        h = 0
    end

    data = zeros(Uint64, length(seq.data))
    j = length(data)
    @inbounds for i in 1:length(data)
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


# Algorithms
# ----------


# compute the mismatch count between two kmers stored in x, y
function nucmismatches(x::Uint64, y::Uint64)
    xyxor = x $ y
    return count_ones((xyxor & 0x5555555555555555) | ((xyxor & 0xAAAAAAAAAAAAAAAA) >>> 1))
end


# return a mask of the first k bits
function makemask(k::Integer)
    return 0xffffffffffffffff >> (64 - k)
end


function mismatches{T}(a::NucleotideSequence{T}, b::NucleotideSequence{T},
                       nmatches::Bool=false)
    # we need to assume that a is aligned with its data
    if a.part.start != 1
        if b.part.start == 1
            return mismatches(b, a)
        else
            if length(a) < b
                orphan!(a)
            else
                orphan!(b)
                return mismatches(b, a)
            end
        end
    end

    d0, r0 = divrem(b.part.start - 1, 64)
    h = 64 - r0
    k = r0

    hmask = makemask(h)
    kmask = makemask(k)
    count = 0
    j = d0 + 1
    for i in 1:length(a.data)
        count += nucmismatches(a.data[i] & hmask, b.data[j] >>> k)
        j += 1
        if j > length(b.data)
            break
        end
        count += nucmismatches(a.data[i] & kmask, b.data[j] << h)
    end

    nsa = ns(a)
    nsb = ns(b)
    nsa_state = start(nsa)
    nsb_state = start(nsb)
    a_done = done(nsa, nsa_state)
    b_done = done(nsb, nsb_state)
    i, nsa_state = a_done ? (0, nsa_state) : next(nsa, nsa_state)
    j, nsb_state = b_done ? (0, nsb_state) : next(nsb, nsb_state)

    while true
        a_hasnext = !done(nsa, nsa_state)
        b_hasnext = !done(nsb, nsb_state)

        if !a_done && (b_done || i < j)
            nucsmatch = getnuc(T, a.data, i + a.part.start - 1) ==
                        getnuc(T, b.data, i + b.part.start - 1)
            if nucsmatch
                if !nmatches
                    count += 1
                end
            else
                if nmatches
                    count -= 1
                end
            end

            if a_hasnext; i, nsa_state = next(nsa, nsa_state)
            else; a_done = true; end
        elseif !b_done && (a_done || j < i)
            nucsmatch = getnuc(T, a.data, j + a.part.start - 1) ==
                        getnuc(T, b.data, j + b.part.start - 1)
            if nucsmatch
                if !nmatches
                    count += 1
                end
            else
                if nmatches
                    count -= 1
                end
            end

            if b_hasnext; j, nsb_state = next(nsb, nsb_state)
            else; b_done = true; end
        elseif !a_done && !b_done
            nucsmatch = getnuc(T, a.data, i + a.part.start - 1) ==
                        getnuc(T, b.data, i + b.part.start - 1)
            if !nucsmatch
                count -= 1
            end

            if a_hasnext; i, nsa_state = next(nsa, nsa_state)
            else; a_done = true; end
            if b_hasnext; j, nsb_state = next(nsb, nsb_state)
            else; b_done = true; end
        else
            break
        end
    end

    return count
end


# Call a function on every evenly spaces subsequence of a given length.
#
# This avoids a little overhead by not reallocating subsequences but rather
# modifying one in place and passing it to the function repeatedly.
#
# Args:
#   f: a function of the form f(subseq::NucleotideSequence)
#   seq: a sequence
#   len: length of the subsequences
#   step: step between the start position of each subsequenc
#
function eachsubseq{T}(f::Base.Callable, seq::NucleotideSequence{T}, len::Integer,
                       step::Integer=1)
    subseq = NucleotideSequence{T}(seq)
    for i in 1:step:(length(seq) - len + 1)
        start = seq.part.start + i - 1
        subseq.part = start:(start + len - 1)
        f(subseq)
    end
end


# Call a function on every k-mer.
#
#
function eachkmer{T}(f::Base.Callable, seq::NucleotideSequence{T}, k::Integer)
    if k > 32
        error("k must be ≤ 32 in eachkmer")
    end
    x = uint64(0)
    mask = makemask(2 * k)
    skip = k - 1
    len = length(seq)
    shift = 2 * (k - 1)
    d, r = divrem(seq.part.start - 1, 64)

    # manually iterate through N positions
    ns_it = ns(seq)
    ns_it_state = start(ns_it)
    if done(ns_it, ns_it_state)
        next_n_pos = length(seq) + 1
    else
        next_n_pos, ns_it_state = next(ns_it, ns_it_state)
    end

    i = 1
    while i <= len
        # skip over any kmer containing an N
        if i == next_n_pos
            if done(ns_it, ns_it_state)
                next_n_pos = length(seq) + 1
            else
                next_n_pos, ns_it_state = next(ns_it, ns_it_state)
            end
            skip = k
        end

        x = (x >>> 2) | (((seq.data[d + 1] >>> r) & 0b11) << shift)

        r += 2
        if r >= 64
            r = 0
            d += 1
        end

        if skip == 0
            f(i, x & mask)
            skip = 1
        else
            skip -= 1
        end
        i += 1
    end
end



# Nucleotide Composition
# ----------------------


# Count A, C, T/U, G respectively in a kmer stored in a Uint64
function count_a(x::Uint64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & 0x5555555555555555)
end

function count_c(x::Uint64)
    return count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
end

function count_g(x::Uint64)
    return count_ones(((x >>> 1) & (~x)) & 0x5555555555555555)
end

function count_t(x::Uint64)
    return count_ones((x & (x >>> 1)) & 0x5555555555555555)
end


type NucleotideCounts{T <: Nucleotide}
    a::Uint
    c::Uint
    g::Uint
    t::Uint
    n::Uint

    function NucleotideCounts()
        new(0, 0, 0, 0, 0)
    end
end


typealias DNANucleotideCounts NucleotideCounts{DNANucleotide}
typealias RNANucleotideCounts NucleotideCounts{RNANucleotide}


function getindex{T}(counts::NucleotideCounts{T}, nt::T)
    return getfield(counts, int(convert(Uint, nt) + 1))
end


function setindex!{T}(counts::NucleotideCounts{T}, d::Integer, nt::T)
    return setfield!(counts, int(convert(Uint, nt) + 1), counts[nt] + 1)
end


# pad strings so they are right justified when printed
function format_counts(xs)
    strings = String[string(x) for x in xs]
    len = maximum(map(length, strings))
    if i in 1:length(strings)
        strings[i] = string(repeat(" ", len - length(strings[i])), strings[i])
    end
    return strings
end


function show(io::IO, counts::DNANucleotideCounts)
    count_strings = format_counts(
        [counts[DNA_A], counts[DNA_C], counts[DNA_G], counts[DNA_T], counts[DNA_n]])

    write(io,
        """
        DNANucleotideCounts:
          A => $(count_strings[1])
          C => $(count_strings[2])
          G => $(count_strings[3])
          T => $(count_strings[4])
          N => $(count_strings[5])
        """)
end


function show(io::IO, counts::RNANucleotideCounts)
    count_strings = format_counts(
        [counts[RNA_A], counts[RNA_C], counts[RNA_G], counts[RNA_U], counts[RNA_n]])

    write(io,
        """
        RNANucleotideCounts:
          A => $(count_strings[1])
          C => $(count_strings[2])
          G => $(count_strings[3])
          U => $(count_strings[4])
          N => $(count_strings[5])
        """)
end


function nucleotide_count{T}(seq::NucleotideSequence{T})
    d, r = divrem(2 * (seq.part.start - 1), 64)
    i = 1

    counts = NucleotideCounts{T}()

    # count leading unaligned bases
    for i in 1:r
        counts[getnuc(T, seq.data, seq.part.start + i - 1)] += 1
        i += 1
    end
    if r > 0
        d += 1
    end

    # count aligned bases
    while i + 32 <= length(seq)
        counts[DNA_A] += count_a(seq.data[d + 1])
        counts[DNA_C] += count_c(seq.data[d + 1])
        counts[DNA_G] += count_g(seq.data[d + 1])
        counts[DNA_T] += count_t(seq.data[d + 1])
        d += 1
        i += 32
    end

    # count trailing unaligned bases
    while i <= length(seq)
        counts[getnuc(T, seq.data, seq.part.start + i - 1)] += 1
        i += 1
    end

    # process Ns
    for i in ns(seq)
        counts[getnuc(T, seq.data, seq.part.start + i - 1)] -= 1
        counts[nnucleotide(T)] += 1
    end

    return counts
end


end # module Seq

