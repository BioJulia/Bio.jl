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
    const maxcount = 50
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
        error(BoundsError())
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
end

# Nucleotide complement
function complement(seq::NucleotideSequence)
    seq = copy(seq)
    unsafe_complement!(seq)
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
function reverse_complement(seq::NucleotideSequence)
    seq = reverse(seq)
    unsafe_complement!(seq)
    return seq
end


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
