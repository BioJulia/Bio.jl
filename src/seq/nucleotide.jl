# Nucleotides
# ===========


# Single nucleotides are represented in bytes using just the two low-order bits
abstract Nucleotide

bitstype 8 DNANucleotide <: Nucleotide
bitstype 8 RNANucleotide <: Nucleotide


# Conversion from/to integers
# ---------------------------

convert(::Type{DNANucleotide}, nt::Uint8) = box(DNANucleotide, unbox(Uint8, nt))
convert(::Type{Uint8}, nt::DNANucleotide) = box(Uint8, unbox(DNANucleotide, nt))

convert(::Type{RNANucleotide}, nt::Uint8) = box(RNANucleotide, unbox(Uint8, nt))
convert(::Type{Uint8}, nt::RNANucleotide) = box(Uint8, unbox(RNANucleotide, nt))

convert{T<:Unsigned, S<:Nucleotide}(::Type{T}, nt::S) = box(T, Base.zext_int(T, unbox(S, nt)))
convert{T<:Unsigned, S<:Nucleotide}(::Type{S}, nt::T) = convert(S, convert(Uint8, nt))


# Nucleotide encoding definition
# ------------------------------

const DNA_A = convert(DNANucleotide, 0b000)
const DNA_C = convert(DNANucleotide, 0b001)
const DNA_G = convert(DNANucleotide, 0b010)
const DNA_T = convert(DNANucleotide, 0b011)
const DNA_N = convert(DNANucleotide, 0b100)
const DNA_INVALID = convert(DNANucleotide, 0b1000) # Indicates invalid DNA when converting string

nnucleotide(::Type{DNANucleotide}) = DNA_N

const RNA_A = convert(RNANucleotide, 0b000)
const RNA_C = convert(RNANucleotide, 0b001)
const RNA_G = convert(RNANucleotide, 0b010)
const RNA_U = convert(RNANucleotide, 0b011)
const RNA_N = convert(RNANucleotide, 0b100)
const RNA_INVALID = convert(RNANucleotide, 0b1000) # Indicates invalid DNA when converting string

nnucleotide(::Type{RNANucleotide}) = RNA_N


# Conversion from Char
# --------------------

# lookup table for characters in 'A':'t'
const char_to_dna = [
     DNA_A,       DNA_INVALID, DNA_C,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_N,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_T,       DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_A,       DNA_INVALID, DNA_C,       DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_G,       DNA_INVALID, DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_N,       DNA_INVALID, DNA_INVALID,
     DNA_INVALID, DNA_INVALID, DNA_INVALID, DNA_T ]


function convert(::Type{DNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 't' ? char_to_dna[c - 'A' + 1] : DNA_INVALID
    @assert nt != DNA_INVALID error("$(c) is not a valid DNA nucleotide")
    return nt
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
    @assert nt != RNA_INVALID error("$(c) is not a valid RNA nucleotide")
    return nt
end


# Conversion to Char
# ------------------

const dna_to_char = ['A', 'C', 'G', 'T', 'N']
convert(::Type{Char}, nt::DNANucleotide) = dna_to_char[convert(Uint8, nt) + 1]

const rna_to_char = ['A', 'C', 'G', 'U', 'N']
convert(::Type{Char}, nt::RNANucleotide) = rna_to_char[convert(Uint8, nt) + 1]


# Basic functions
# ---------------

function show{T<:Nucleotide}(io::IO, nt::T)
    if nt == DNA_INVALID
        write(io, "Invalid DNA Nucleotide")
    elseif nt == RNA_INVALID
        write(io, "Invalid RNA Nucleotide")
    else
        write(io, convert(Char, nt))
    end
end



# Nucleotide Sequence
# ===================

# A general purpose nucleotide sequence representation
#
# Nucleotide sequences are 2-bit encoded and packed into Uint64s. 'N's are
# represented with an N mask BitArray. If the ns[i] bit is set, then the
# sequence may have any nucleotide at that position and it must be ignored.
#
# Sequences are immutable by convention. Low-level functions will mutate the
# underlying data, any user-facing function that does so should have `unsafe` in
# the name. Immutability allows two sequences to share the same underlying data,
# so that subsequences can be initialized without significant copying or
# allocation.

type NucleotideSequence{T<:Nucleotide}
    data::Vector{Uint64} # 2-bit encoded sequence
    ns::BitVector        # 'N' mask
    part::UnitRange{Int} # interval within `data` and `ns` defining the (sub)sequence
end


# Constructors
# ------------

# Construct an empty sequence
NucleotideSequence{T<:Nucleotide}(::Type{T}) = NucleotideSequence{T}(zeros(Uint64, 0), BitVector(0), 1:0)

# Construct a subsequence of another nucleotide sequence
function NucleotideSequence{T<:Nucleotide}(::Type{T}, other::NucleotideSequence, part::UnitRange)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    if start < other.part.start || stop > other.part.stop
        error("Invalid subsequence range")
    end
    return NucleotideSequence{T}(other.data, other.ns, part)
end

# Return the number of Uint64s needed to represent a sequence of length n
function seq_data_len(n::Integer)
    d, r = divrem(n, 32)
    return d + (r > 0 ? 1 : 0)
end

# Construct a sequence from a string
function NucleotideSequence{T<:Nucleotide}(::Type{T}, seq::String)
    len  = seq_data_len(length(seq))
    data = zeros(Uint64, len)
    ns   = BitArray(length(seq))
    fill!(ns, false)

    j = start(seq)
    idx = 1
    @inbounds for i in 1:length(data)
        shift = 0
        while shift < 64 && !done(seq, j)
            c, j = next(seq, j)
            nt   = convert(T, c)
            if nt == nnucleotide(T)
                ns[idx] = true
            else
                data[i] |= convert(Uint64, nt) << shift
            end
            idx += 1
            shift += 2
        end
    end

    return NucleotideSequence{T}(data, ns, 1:length(seq))
end

# Aliases and contructors
# -----------------------

typealias DNASequence NucleotideSequence{DNANucleotide}
DNASequence() = NucleotideSequence(DNANucleotide)
DNASequence(other::NucleotideSequence, part::UnitRange) = NucleotideSequence(DNANucleotide, other, part)
DNASequence(seq::String) = NucleotideSequence(DNANucleotide, seq)

typealias RNASequence NucleotideSequence{RNANucleotide}
RNASequence() = NucleotideSequence(RNANucleotide)
RNASequence(other::NucleotideSequence, part::UnitRange) = NucleotideSequence(RNANucleotide, other, part)
RNASequence(seq::String) = NucleotideSequence(RNANucleotide, seq)


# Conversion
# ----------

# Convert from/to Strings
convert(::Type{DNASequence}, seq::String) = DNASequence(seq)
convert(::Type{RNASequence}, seq::String) = RNASequence(seq)
convert(::Type{String}, seq::NucleotideSequence) = convert(String, [convert(Char, x) for x in seq])


# Convert between RNA and DNA
convert(::Type{RNASequence}, seq::DNASequence) = RNASequence(seq.data, seq.ns, seq.part)
convert(::Type{DNASequence}, seq::RNASequence) = DNASequence(seq.data, seq.ns, seq.part)


# Basic Functions
# ---------------

length(seq::NucleotideSequence) = length(seq.part)
endof(seq::NucleotideSequence)  = length(seq)

# Pretting printing of sequences.
function show{T}(io::IO, seq::NucleotideSequence{T})
    len = length(seq)
    write(io, "$(string(len))nt ", T == DNANucleotide ? "DNA" : "RNA", " Sequence\n ")

    # don't show more than this many characters to avoid filling the screen
    # with junk
    const maxcount = Base.tty_size()[2] - 2
    if len > maxcount
        for nt in seq[1:div(maxcount, 2) - 1]
            write(io, convert(Char, nt))
        end
        write(io, "…")
        for nt in seq[(end - (div(maxcount, 2) - 1)):end]
            write(io, convert(Char, nt))
        end
    else
        for nt in seq
            write(io, convert(Char, nt))
        end
    end
end

function =={T}(a::NucleotideSequence{T}, b::NucleotideSequence{T})
    if a.data === b.data && a.ns === b.ns && a.part == b.part
        return true
    end

    if length(a) != length(b)
        return false
    end

    for (u, v) in zip(a, b)
        if u != v
            return false
        end
    end

    return true
end

# Get the nucleotide at position i, ignoring the N mask.
function getnuc(T::Type, data::Vector{Uint64}, i::Integer)
    d, r = divrem(i - 1, 32)
    return convert(T, convert(Uint8, (data[d + 1] >>> (2*r)) & 0b11))
end

function getindex{T}(seq::NucleotideSequence{T}, i::Integer)
    if i > length(seq) || i < 1
        throw(BoundsError())
    end
    i += seq.part.start - 1
    if seq.ns[i]
        return nnucleotide(T)
    else
        return getnuc(T, seq.data, i)
    end
end

# Construct a subesequence
getindex{T}(seq::NucleotideSequence{T}, r::UnitRange) = NucleotideSequence{T}(seq, r)


# Replace a NucleotideSequence's data with a copy, copying only what's needed.
#
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

    data   = zeros(Uint64, seq_data_len(length(seq)))
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
    seq.ns   = seq.ns[seq.part.start:seq.part.stop]
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
copy{T}(seq::NucleotideSequence{T}) = orphan!(NucleotideSequence{T}(seq.data, seq.ns, seq.part), true)


# Iterating throug nucleotide sequences
# TODO: This can be made a lot faster.
start(seq::NucleotideSequence) = seq.part.start - 1

function next{T}(seq::NucleotideSequence{T}, i)
    # Check bounds
    (seq.part.start - 1 <= i < seq.part.stop) || throw(BoundsError())

    nvalue, _ = next(seq.ns, i)
    if nvalue
        return (nnucleotide(T), i + 1)
    else
        return (getnuc(T, seq.data, i + 1), i + 1)
    end
end

done(seq::NucleotideSequence, i) = i >= seq.part.stop

# String Decorator
# ----------------

# Enable building sequence literals like: dna"ACGTACGT" and rna"ACGUACGU"
macro dna_str(seq, flags...)
    return DNASequence(seq)
end

macro rna_str(seq, flags...)
    return RNASequence(seq)
end


# Transformations
# ---------------

# In-place complement (nucleotide complement is equivalent to bitwise complement
# in the encoding used)
function unsafe_complement!(seq::NucleotideSequence)
    @inbounds for i in 1:length(seq.data)
        seq.data[i] = ~seq.data[i]
    end
    return seq
end

# Nucleotide complement
complement(seq::NucleotideSequence) = unsafe_complement!(copy(seq))

# Nucleotide reverse. Reverse a kmer stored in a Uint64.
function nucrev(x::Uint64)
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>>  2
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >>> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >>> 32
     return x
end

# Return a reversed copy of seq
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

# Return the reverse complement of seq
reverse_complement(seq::NucleotideSequence) = unsafe_complement!(reverse(seq))


# SequenceNIterator
# -----------------

# Iterate through positions in the sequence with Ns
#
# This can be much faster than testing every position (seq.ns[i]) since
# it can skip over 64 positions at time if they don't have 'N's.
immutable SequenceNIterator
    ns::BitVector
    part::UnitRange{Int}
end

SequenceNIterator(seq::NucleotideSequence) = SequenceNIterator(seq.ns, seq.part)
npositions(seq::NucleotideSequence)        = SequenceNIterator(seq)


# Find the next N in the sequence starting at position i.
#
# Return any position past the end of the sequence if there are no more Ns.
#
function nextn(it::SequenceNIterator, i)
    d, r = divrem(i - 1, 64)
    while d < length(it.ns.chunks) && it.ns.chunks[d + 1] >>> r == 0 && d * 64 < it.part.stop
        d += 1
        r  = 0
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


# Iterating through SequenceNIterator
# -----------------------------------

start(it::SequenceNIterator) = nextn(it, it.part.start)

function next(it::SequenceNIterator, i)
    d, r = divrem(i - 1, 64)
    next_i = nextn(it, i + 1)
    return i + it.part.start - 1, next_i
end

done(it::SequenceNIterator, i) = i > it.part.stop

function hasn(seq::NucleotideSequence)
    it = npositions(seq)
    return !done(it, start(it))
end

# TODO: Implement length for SequenceNIterators to use comprehensions


# Mismatch counting
# -----------------

# compute the mismatch count between two kmers stored in x, y
function nucmismatches(x::Uint64, y::Uint64)
    xyxor = x $ y
    return count_ones((xyxor & 0x5555555555555555) | ((xyxor & 0xAAAAAAAAAAAAAAAA) >>> 1))
end

# return a mask of the first k bits of a Uint64
function makemask(k::Integer)
    return 0xffffffffffffffff >> (64 - k)
end

# Return the number of mismatches between a and b.
#
# If a and b are of differing lengths, only the first min(length(a), length(b))
# nucleotides are compared.
#
# Args:
#   a: first sequence to compare
#   b: second sequence to compare
#   nbatches: if true, N matches anything, if false, N matches only itself.
#
# Returns:
#   The number of mismatches.
#
function mismatches{T}(a::NucleotideSequence{T}, b::NucleotideSequence{T},
                       nmatches::Bool=false)

    # we need to assume that a is aligned with its data, rearrange the
    # comparison if that's not the case.
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

    # Count mismatches, ignoring the presence of 'N's in the sequence.
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

    # Here's the ugly part. We've just counted mismtaches without taking the N
    # mask into account. If 'N's are present in the sequence, that mismatch
    # count may be too high or too low, so we walk through all N positions in
    # both sequences in unision, correcting the count where needed.
    nsa       = npositions(a)
    nsb       = npositions(b)
    nsa_state = start(nsa)
    nsb_state = start(nsb)
    a_done    = done(nsa, nsa_state)
    b_done    = done(nsb, nsb_state)

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



# K-mer
# =====


# A Kmer is a sequence <= 32nt, without any 'N's, packed in a single 64 bit value.
#
# While NucleotideSequence is an efficient general-purpose sequence
# representation, Kmer is useful for applications like assembly, k-mer counting,
# k-mer based quantification in RNA-Seq, etc that rely on manipulating many
# short sequences as efficiently (space and time) as possible.

bitstype 64 Kmer{T<:Nucleotide, K}

typealias DNAKmer{K} Kmer{DNANucleotide, K}
typealias RNAKmer{K} Kmer{RNANucleotide, K}
typealias Codon RNAKmer{3}


# Conversion
# ----------

# Conversion to/from Uint64
convert{K}(::Type{DNAKmer{K}}, x::Uint64) = box(DNAKmer{K}, unbox(Uint64, x))
convert{K}(::Type{RNAKmer{K}}, x::Uint64) = box(RNAKmer{K}, unbox(Uint64, x))
convert{K}(::Type{Uint64}, x::DNAKmer{K})       = box(Uint64, unbox(DNAKmer{K}, x))
convert{K}(::Type{Uint64}, x::RNAKmer{K})       = box(Uint64, unbox(RNAKmer{K}, x))


# Convert from String
function convert{T, K}(::Type{Kmer{T, K}}, seq::String)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a string of length $(length(seq))")

    x     = @compat UInt64(0)
    shift = 0
    for (i, c) in enumerate(seq)
        nt = convert(T, c)
        @assert nt != nnucleotide(T) error("A K-mer may not contain an N in its sequence")

        x |= convert(Uint64, nt) << shift
        shift += 2
    end

    return convert(Kmer{T, K}, x)
end

convert{T}(::Type{Kmer{T}}, seq::String) = convert(Kmer{T, length(seq)}, seq)

# Convert to String
convert{T, K}(::Type{String}, seq::Kmer{T, K}) = convert(String, [convert(Char, x) for x in seq])


# Convert from NucleotideSequence
function convert{T, K}(::Type{Kmer{T, K}}, seq::NucleotideSequence{T})
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    @assert length(seq) == K error("Cannot construct a $(K)-mer from a NucleotideSequence of length $(length(seq))")

    x     = @compat UInt64(0)
    shift = 0
    for (i, nt) in enumerate(seq)
        if nt == nnucleotide(T)
            error("A K-mer may not contain an N in its sequence")
        end
        x |= convert(Uint64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end

convert{T}(::Type{Kmer}, seq::NucleotideSequence{T})    = convert(Kmer{T, length(seq)}, seq)
convert{T}(::Type{Kmer{T}}, seq::NucleotideSequence{T}) = convert(Kmer{T, length(seq)}, seq)

# Convert to NucleotideSequence
function convert{T, K}(::Type{NucleotideSequence{T}}, x::Kmer{T, K})
    ns = BitVector(K)
    fill!(ns, false)
    return NucleotideSequence{T}([convert(Uint64, x)], ns, 1:K)
end

convert{T, K}(::Type{NucleotideSequence}, x::Kmer{T, K}) = convert(NucleotideSequence{T}, x)




# Constructors
# ------------

# From strings
dnakmer(seq::String) = convert(DNAKmer, seq)
rnakmer(seq::String) = convert(RNAKmer, seq)

# Constructors taking a sequence of nucleotides
function kmer{T <: Nucleotide}(nts::T...)
    K = length(nts)
    if K > 32
        error(string("Cannot construct a K-mer longer than 32nt."))
    end

    x = @compat UInt64(0)
    shift = 0
    for (i, nt) in enumerate(nts)
        if nt == nnucleotide(T)
            error("A Kmer may not contain an N in its sequence")
        end
        x |= convert(Uint64, nt) << shift
        shift += 2
    end
    return convert(Kmer{T, K}, x)
end

function kmer(seq::DNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(DNAKmer{length(seq)}, seq)
end

function kmer(seq::RNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(RNAKmer{length(seq)}, seq)
end

# call kmer with @inline macro would reduce the performance significantly?
function dnakmer(seq::DNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(DNAKmer{length(seq)}, seq)
end

function rnakmer(seq::RNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(RNAKmer{length(seq)}, seq)
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

function getindex{T, K}(x::Kmer{T, K}, i::Integer)
    if i < 1 || i > K
        throw(BoundsError())
    end
    return convert(T, (convert(Uint64, x) >>> (2*(i-1))) & 0b11)
end


function show{K}(io::IO, x::DNAKmer{K})
    write(io, "DNA $(K)-mer:\n ")
    for i in 1:K
        write(io, convert(Char, x[i]))
    end
end


function show{K}(io::IO, x::RNAKmer{K})
    write(io, "RNA $(K)-mer:\n ")
    for i in 1:K
        write(io, convert(Char, x[i]))
    end
end


isless{T, K}(x::Kmer{T, K}, y::Kmer{T, K}) = convert(Uint64, x) < convert(Uint64, y)

length{T, K}(x::Kmer{T, K}) = K

# Iterating over nucleotides
start(x::Kmer) = 1

function next{T, K}(x::Kmer{T, K}, i::Int)
    nt = convert(T, (convert(Uint64, x) >>> (2*(i-1))) & 0b11)
    return (nt, i + 1)
end

done{T, K}(x::Kmer{T, K}, i::Int) = i > K


# Other functions
# ---------------

reverse{T, K}(x::Kmer{T, K}) = convert(Kmer{T, K}, nucrev(convert(Uint64, x)) >>> (2 * (32 - K)))

function complement{T, K}(x::Kmer{T, K})
    return convert(Kmer{T, K},
        (~convert(Uint64, x)) & (0xffffffffffffffff >>> (2 * (32 - K))))
end

reverse_complement{T, K}(x::Kmer{T, K}) = complement(reverse(x))

mismatches{T, K}(x::Kmer{T, K}, y::Kmer{T, K}) = nucmismatches(convert(Uint64, x), convert(Uint64, y))


# A canonical k-mer is the numerical lesser of a k-mer and its reverse complement.
# This is useful in hashing/counting k-mers in data that is not strand specific,
# and thus observing k-mer is equivalent to observing its reverse complement.
function canonical{T, K}(x::Kmer{T, K})
    y = reverse_complement(x)
    return x < y ? x : y
end




# EachKmerIterator and EachKmerIteratorState
# ==========================================

# Iterate through every k-mer in a nucleotide sequence
immutable EachKmerIterator{T, K}
    seq::NucleotideSequence{T}
    nit::SequenceNIterator
    step::Int
end


immutable EachKmerIteratorState{T, K}
    i::Int
    x::Uint64
    next_n_pos::Int
    nit_state::Int
end


function eachkmer{T}(seq::NucleotideSequence{T}, k::Integer, step::Integer=1)
    if k < 0
        error("K must be ≥ 0 in EachKmer")
    elseif k > 32
        error("K must be ≤ 32 in EachKmer")
    end

    if step < 1
        error("step must be ≥ 1")
    end

    return EachKmerIterator{T, k}(seq, npositions(seq), step)
end


function nextkmer{T, K}(it::EachKmerIterator{T, K},
                        state::EachKmerIteratorState{T, K}, skip::Int)
    i = state.i + 1
    x = state.x
    next_n_pos = state.next_n_pos
    nit_state = state.nit_state

    shift = 2 * (K - 1)
    d, r = divrem(2 * (it.seq.part.start + i - 2), 64)
    while i <= length(it.seq)
        while next_n_pos < i
            if done(it.nit, nit_state)
                next_n_pos = length(it.seq) + 1
                break
            else
                next_n_pos, nit_state = next(it.nit, nit_state)
            end
        end

        if i - K + 1 <= next_n_pos <= i
            off = it.step * @compat ceil(Int, (K - skip) / it.step)
            if skip < K
                skip += off
            end
        end

        x = (x >>> 2) | (((it.seq.data[d + 1] >>> r) & 0b11) << shift)

        if skip == 0
            break
        end
        skip -= 1

        r += 2
        if r == 64
            r = 0
            d += 1
        end
        i += 1
    end

    return EachKmerIteratorState{T, K}(i, x, next_n_pos, nit_state)
end


function start{T, K}(it::EachKmerIterator{T, K})
    nit_state = start(it.nit)
    if done(it.nit, nit_state)
        next_n_pos = length(it.seq) + 1
    else
        next_n_pos, nit_state = next(it.nit, nit_state)
    end

    state = EachKmerIteratorState{T, K}(0, (@compat UInt64(0)), next_n_pos, nit_state)
    return nextkmer(it, state, K - 1)
end


function next{T, K}(it::EachKmerIterator{T, K},
                    state::EachKmerIteratorState{T, K})
    value = convert(Kmer{T, K}, state.x)
    next_state = nextkmer(it, state, it.step - 1)
    return (state.i - K + 1, value), next_state
end


function done{T, K}(it::EachKmerIterator{T, K},
                    state::EachKmerIteratorState{T, K})
    return state.i > length(it.seq)
end

# TODO: count_nucleotides



# Nucleotide Composition
# ======================

type NucleotideCounts{T <: Nucleotide}
    a::Uint
    c::Uint
    g::Uint
    t::Uint # also hold 'U' count when T == RNANucleotide
    n::Uint

    function NucleotideCounts()
        new(0, 0, 0, 0, 0)
    end
end


# Aliases
# -------

typealias DNANucleotideCounts NucleotideCounts{DNANucleotide}
typealias RNANucleotideCounts NucleotideCounts{RNANucleotide}


# Constructors
# ------------

# Count A, C, T/U, G respectively in a kmer stored in a Uint64
function count_a(x::Uint64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & 0x5555555555555555)
end
count_c(x::Uint64) = count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
count_g(x::Uint64) = count_ones(((x >>> 1) & (~x)) & 0x5555555555555555)
count_t(x::Uint64) = count_ones((x    & (x >>> 1)) & 0x5555555555555555)

function NucleotideCounts{T}(seq::NucleotideSequence{T})
    dn, rn = divrem(seq.part.start - 1, 64)

    d = 2*dn
    r = 2*dn

    i = 1
    counts = NucleotideCounts{T}()

    # count leading unaligned bases
    for i in 1:r
        counts[seq[i]] += 1
        i += 1
    end
    if r > 0
        d += 1
    end

    # maybe just skip over blocks of Ns as I go?
    while i + 63 <= length(seq)
        # handle the all-N case
        if seq.ns.chunks[dn + 1] == 0xffffffffffffffff
            counts.n += 64
        else
            counts.a += count_a(seq.data[d + 1]) + count_a(seq.data[d + 2])
            counts.c += count_c(seq.data[d + 1]) + count_c(seq.data[d + 2])
            counts.g += count_g(seq.data[d + 1]) + count_g(seq.data[d + 2])
            counts.t += count_t(seq.data[d + 1]) + count_t(seq.data[d + 2])

            x = seq.ns.chunks[dn + 1]
            if x != 0
                for j in 1:64
                    if x & 0x01 != 0
                        counts.n += 1
                        counts[getnuc(T, seq.data, seq.part.start + i + j - 2)] -= 1
                    end

                    x >>= 1
                    if x == 0
                        break
                    end
                end
            end
        end

        dn += 1
        d += 2
        i += 64
    end

    # count trailing unaligned bases
    while i <= length(seq)
        counts[seq[i]] += 1
        i += 1
    end

    return counts
end

# Construct from K-mers
function NucleotideCounts{T,K}(seq::Kmer{T, K})
    x         = convert(Uint64, seq)
    counts    = NucleotideCounts{T}()
    counts.a += count_a(x) - 32 + K # Take leading zeros into account
    counts.c += count_c(x)
    counts.g += count_g(x)
    counts.t += count_t(x)
    return counts
end

# Basic Functions
# ---------------

getindex{T}(counts::NucleotideCounts{T}, nt::T) = getfield(counts, (@compat Int(convert(Uint, nt) + 1)))
setindex!{T}(counts::NucleotideCounts{T}, c::Integer, nt::T) = setfield!(counts, (@compat Int(convert(Uint, nt) + 1)), c)

# Pad strings so they are right justified when printed
function format_counts(xs)
    strings = String[string(x) for x in xs]
    len = maximum(map(length, strings))
    for i in 1:length(strings)
        strings[i] = string(repeat(" ", len - length(strings[i])), strings[i])
    end
    return strings
end


# Pretty printing of NucleotideCounts
function show(io::IO, counts::DNANucleotideCounts)
    count_strings = format_counts(
        [counts[DNA_A], counts[DNA_C], counts[DNA_G], counts[DNA_T], counts[DNA_N]])

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
        [counts[RNA_A], counts[RNA_C], counts[RNA_G], counts[RNA_U], counts[RNA_N]])

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
