# Nucleotides
# ===========


# Single nucleotides are represented in bytes using just the two low-order bits
abstract Nucleotide
bitstype 8 DNANucleotide <: Nucleotide
bitstype 8 RNANucleotide <: Nucleotide


# Conversion from/to integers
# ---------------------------

"Convert a Uint8 to a DNANucleotide"
convert(::Type{DNANucleotide}, nt::Uint8) = box(DNANucleotide, unbox(Uint8, nt))

"Convert a DNANucleotide to a Uint8"
convert(::Type{Uint8}, nt::DNANucleotide) = box(Uint8, unbox(DNANucleotide, nt))

"Convert a Uint8 to a RNANucleotide"
convert(::Type{RNANucleotide}, nt::Uint8) = box(RNANucleotide, unbox(Uint8, nt))

"Convert a RNANucleotide to a Uint8"
convert(::Type{Uint8}, nt::RNANucleotide) = box(Uint8, unbox(RNANucleotide, nt))

"Convert a RNANucleotide to a Uint8"
convert{T<:Unsigned, S<:Nucleotide}(::Type{T}, nt::S) = box(T, Base.zext_int(T, unbox(S, nt)))

"Convert a RNANucleotide to a Uint8"
convert{T<:Unsigned, S<:Nucleotide}(::Type{S}, nt::T) = convert(S, convert(Uint8, nt))

# Nucleotide encoding definition
# ------------------------------

# DNA Nucleotides

"DNA Adenine"
const DNA_A = convert(DNANucleotide, 0b000)

"DNA Cytosine"
const DNA_C = convert(DNANucleotide, 0b001)

"DNA Guanine"
const DNA_G = convert(DNANucleotide, 0b010)

"DNA Thymine"
const DNA_T = convert(DNANucleotide, 0b011)

"DNA Any Nucleotide"
const DNA_N = convert(DNANucleotide, 0b100)

"DNA Invalid Nucleotide"
const DNA_INVALID = convert(DNANucleotide, 0b1000) # Indicates invalid DNA when converting string

"Returns Any DNA Nucleotide (DNA_N)"
nnucleotide(::Type{DNANucleotide}) = DNA_N

# RNA Nucleotides

"RNA Adenine"
const RNA_A = convert(RNANucleotide, 0b000)

"RNA Cytosine"
const RNA_C = convert(RNANucleotide, 0b001)

"RNA Guanine"
const RNA_G = convert(RNANucleotide, 0b010)

"RNA Uracil"
const RNA_U = convert(RNANucleotide, 0b011)

"Any RNA Nucleotide"
const RNA_N = convert(RNANucleotide, 0b100)

"Invalid RNA Nucleotide"
const RNA_INVALID = convert(RNANucleotide, 0b1000) # Indicates invalid RNA when converting string

"Returns Any RNA Nucleotide (RNA_N)"
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

"Convert a Char to a DNANucleotide"
function convert(::Type{DNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 't' ? char_to_dna[c - 'A' + 1] : DNA_INVALID
    @assert nt != DNA_INVALID error(" $(c) is not a valid DNA nucleotide")
    return nt
end

function unsafe_ascii_byte_to_nucleotide(T::Type{DNANucleotide}, c::Uint8)
    @inbounds nt = char_to_dna[c - 0x41 + 1]
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

"Convert a Char to a RNANucleotide"
function convert(::Type{RNANucleotide}, c::Char)
    @inbounds nt = 'A' <= c <= 'u' ? char_to_rna[c - 'A' + 1] : RNA_INVALID
    @assert nt != RNA_INVALID error(" $(c) is not a valid RNA nucleotide")
    return nt
end

function unsafe_ascii_byte_to_nucleotide(T::Type{RNANucleotide}, c::Uint8)
    @inbounds nt = char_to_rna[c - 0x41 + 1]
    return nt
end


# Conversion to Char
# ------------------

const dna_to_char = ['A', 'C', 'G', 'T', 'N']

"Convert a DNANucleotide to a Char"
convert(::Type{Char}, nt::DNANucleotide) = dna_to_char[convert(Uint8, nt) + 1]

const rna_to_char = ['A', 'C', 'G', 'U', 'N']

"Convert a RNANucleotide to a Char"
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

"""
`NucleotideSequence(DNANucleotide|RNANucleotide)`

Construct an empty nucleotide sequence of the given type
"""
NucleotideSequence{T<:Nucleotide}(::Type{T}) = NucleotideSequence{T}(zeros(Uint64, 0), BitVector(0), 1:0)


"""
`NucleotideSequence(DNANucleotide|RNANucleotide, other::NucleotideSequence, part::UnitRange)`

Construct a subsequence of the given type from another nucleotide sequence
"""
function NucleotideSequence{T<:Nucleotide}(::Type{T}, other::NucleotideSequence, part::UnitRange)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    if start < other.part.start || stop > other.part.stop
        error("Invalid subsequence range")
    end
    return NucleotideSequence{T}(other.data, other.ns, part)
end


# Faster divrem(n, 32)
function divrem32(n::Integer)
    return (n >> 5, n & 0b11111)
end


# Faster divrem(n, 64)
function divrem64(n::Integer)
    return (n >> 6, n & 0b111111)
end


# Return the number of Uint64s needed to represent a sequence of length n
function seq_data_len(n::Integer)
    d, r = divrem32(n)
    return d + (r > 0 ? 1 : 0)
end


# This is the body of the NucleotideSequence constructor below. It's separated
# into a macro so we can generate two versions depending on wether the `unsafe`
# flag is set.
macro encode_seq(nt_convert_expr)
    quote
        @inbounds begin
            for i in 1:length(data)
                shift = 0
                while shift < 64 && j <= stoppos
                    c = seq[j]
                    j += 1
                    nt = $(nt_convert_expr)
                    if nt == nnucleotide(T)
                        # manually inlined: ns[i] = true
                        d = (idx - 1) >>> 6
                        r = (idx - 1) & 63
                        ns.chunks[d + 1] |= (@compat UInt64(1)) << r
                    else
                        data[i] |= convert(Uint64, nt) << shift
                    end

                    idx += 1
                    shift += 2
                end
            end
        end
    end
end


"""
`NucleotideSequence(DNANucleotide|RNANucleotide, seq::String)`

Construct a subsequence from the `seq` string
"""
function NucleotideSequence{T<:Nucleotide}(::Type{T}, seq::Union(String, Vector{Uint8}),
                                           startpos::Int, stoppos::Int, unsafe::Bool=false)
    len = seq_data_len(stoppos - startpos + 1)
    data = zeros(Uint64, len)
    ns = BitArray(stoppos - startpos + 1)
    fill!(ns, false)

    j = startpos
    idx = 1
    if unsafe
        @encode_seq unsafe_ascii_byte_to_nucleotide(T, c)
    else
        @encode_seq convert(T, convert(Char, c))
    end

    return NucleotideSequence{T}(data, ns, 1:(stoppos - startpos + 1))
end


function NucleotideSequence{T<:Nucleotide}(t::Type{T}, seq::Union(String, Vector{Uint8}))
    return NucleotideSequence(t, seq, 1, length(seq))
end


"""
`NucleotideSequence(chunks::NucleotideSequence...)`

Construct a nucleotide sequence by concatenating the given sequences
"""
function NucleotideSequence{T<:Nucleotide}(chunks::NucleotideSequence{T}...)
    seqlen = 0
    for chunk in chunks
        seqlen += length(chunk)
    end

    datalen = seq_data_len(seqlen)
    data = zeros(Uint64, datalen)
    ns   = BitArray(seqlen)
    newseq = NucleotideSequence{T}(data, ns, 1:seqlen)
    fill!(ns, false)

    pos = 1
    for chunk in chunks
        unsafe_copy!(newseq, pos, chunk)
        pos += length(chunk)
    end

    return newseq
end


(*){T}(chunk1::NucleotideSequence{T}, chunks::NucleotideSequence{T}...) = NucleotideSequence(chunk1, chunks...)


"""
Construct a NucleotideSequence from an array of nucleotides.
"""
function NucleotideSequence{T<:Nucleotide}(seq::AbstractVector{T})
    len = seq_data_len(length(seq))
    data = zeros(Uint64, len)
    ns = BitArray(length(seq))
    fill!(ns, false)

    for (i, nt) in enumerate(seq)
        if nt == nnucleotide(T)
            ns[i] = true
        else
            d, r = divrem32(i - 1)
            data[d + 1] |= convert(Uint64, nt) << (2*r)
        end
    end

    return NucleotideSequence{T}(data, ns, 1:length(seq))
end


"""
`repeat(chunk::NucleotideSequence, n)`

Construct a nucleotide sequence by repeating another sequence `n` times
"""
function repeat{T<:Nucleotide}(chunk::NucleotideSequence{T}, n::Integer)
    seqlen = n * length(chunk)

    datalen = seq_data_len(seqlen)
    data = zeros(Uint64, datalen)
    ns   = BitArray(seqlen)
    newseq = NucleotideSequence{T}(data, ns, 1:seqlen)
    fill!(ns, false)

    pos = 1
    for i in 1:n
        unsafe_copy!(newseq, pos, chunk)
        pos += length(chunk)
    end

    return newseq
end

"Repeat nucleotide sequences"
(^){T}(chunk::NucleotideSequence{T}, n::Integer) = repeat(chunk, n::Integer)


"""
Copy `src` to `dest` starting at position `pos`.

This is unsafe in the following ways:
- Disregards immutability of `dest`
- May write a few bases past `dest[pos + length(src) - 1]`
- Doesn't bounds check anything.

It's really only suitable for use in the concatenation constructor.
"""
function unsafe_copy!{T}(dest::NucleotideSequence{T}, pos::Int, src::NucleotideSequence{T})
    abspos = dest.part.start + pos - 1
    copy!(dest.ns, abspos, src.ns, src.part.start, length(src))

    d1, r1 = divrem(abspos - 1, 32)
    d2, r2 = divrem(src.part.start - 1, 32)

    l = 0
    while l < length(src)
        if r1 == r2 == 0
            dest.data[d1+1] = src.data[d2+1]
        else
            dest.data[d1+1] |= (src.data[d2+1] >> (2*r2)) << (2*r1)
        end

        if r1 > r2
            k = 32 - r1
            d1 += 1
            r1 = 0
            r2 += k
            l += k
        elseif r2 > r1
            k = 32 - r2
            d2 += 1
            r2 = 0
            r1 += k
            l += k
        else # r1 == r2
            k = 32 - r1
            r1 += k
            r2 += k
            if r1 >= 32
                r1 = 0
                r2 = 0
                d1 += 1
                d2 += 1
            end
            l += k
        end
    end

    # zero positions that we've overwritten at the end
    if l > length(src)
        if r1 == 0
            r1 = 32
            d1 -= 1
        end
        dest.data[d1+1] &= 0xffffffffffffffff >> (2 * (32 - r1 + l - length(src)))
    end
end


# Aliases and contructors
# -----------------------

# DNA Sequences
typealias DNASequence NucleotideSequence{DNANucleotide}

"Construct an empty DNA nucleotide sequence"
DNASequence() = NucleotideSequence(DNANucleotide)

"Construct a DNA nucleotide subsequence from another sequence"
DNASequence(other::NucleotideSequence, part::UnitRange) = NucleotideSequence(DNANucleotide, other, part)

"Construct a DNA nucleotide sequence from a String"
DNASequence(seq::String) = NucleotideSequence(DNANucleotide, seq)

"Construct a DNA nucleotide sequence from other sequences"
DNASequence(chunk1::DNASequence, chunks::DNASequence...) = NucleotideSequence(chunk1, chunks...)
DNASequence(seq::Union(Vector{Uint8}, String)) = NucleotideSequence(DNANucleotide, seq)
DNASequence(seq::Union(Vector{Uint8}, String), startpos::Int, endpos::Int, unsafe::Bool=false) = NucleotideSequence(DNANucleotide, seq, startpos, endpos, unsafe)
DNASequence(seq::AbstractVector{DNANucleotide}) = NucleotideSequence(seq)


# RNA Sequences
typealias RNASequence NucleotideSequence{RNANucleotide}

"Construct an empty RNA nucleotide sequence"
RNASequence() = NucleotideSequence(RNANucleotide)

"Construct a RNA nucleotide subsequence from another sequence"
RNASequence(other::NucleotideSequence, part::UnitRange) = NucleotideSequence(RNANucleotide, other, part)

"Construct a RNA nucleotide sequence from a String"
RNASequence(seq::String) = NucleotideSequence(RNANucleotide, seq)

"Construct a RNA nucleotide sequence from other sequences"
RNASequence(chunk1::RNASequence, chunks::RNASequence...) = NucleotideSequence(chunk1, chunks...)
RNASequence(seq::Union(Vector{Uint8}, String)) = NucleotideSequence(RNANucleotide, seq)
RNASequence(seq::Union(Vector{Uint8}, String), startpos::Int, endpos::Int, unsafe::Bool=false) = NucleotideSequence(RNANucleotide, seq, startpos, endpos, unsafe)
RNASequence(seq::AbstractVector{RNANucleotide}) = NucleotideSequence(seq)


# Conversion
# ----------

# Convert from/to Strings

"Convert a String to a DNASequence"
convert(::Type{DNASequence}, seq::String) = DNASequence(seq)

"Convert a String to a RNASequence"
convert(::Type{RNASequence}, seq::String) = RNASequence(seq)

"Convert a NucleotideSequence to a String"
convert(::Type{String}, seq::NucleotideSequence) = convert(String, [convert(Char, x) for x in seq])


# Convert between RNA and DNA

"Convert a DNASequence to a RNASequence"
convert(::Type{RNASequence}, seq::DNASequence) = RNASequence(seq.data, seq.ns, seq.part)

"Convert a RNASequence to a DNASequence"
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
@inline function getnuc(T::Type, data::Vector{Uint64}, i::Integer)
    d, r = divrem32(i - 1)
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
    d0, r0 = divrem32(seq.part.start - 1)

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
@inline function start(seq::NucleotideSequence)
    npos = nextone(seq.ns, seq.part.start)
    return seq.part.start, npos
end

@inline function next{T}(seq::NucleotideSequence{T}, state)
    i, npos = state
    if i == npos
        npos = nextone(seq.ns, i + 1)
        value = nnucleotide(T)
    else
        d, r = divrem32(i - 1)
        @inbounds value = convert(T, ((seq.data[d + 1] >>> (2 * r)) & 0b11) % Uint8)
    end

    return value, (i + 1, npos)
end

@inline done(seq::NucleotideSequence, state) = state[1] > seq.part.stop

# String Decorator
# ----------------

# Enable building sequence literals like: dna"ACGTACGT" and rna"ACGUACGU"
macro dna_str(seq, flags...)
    return DNASequence(remove_newlines(seq))
end

macro rna_str(seq, flags...)
    return RNASequence(remove_newlines(seq))
end

function remove_newlines(seq)
    return replace(seq, r"\r|\n", "")
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


"""
`complement(seq::NucleotideSequence)`

The nucleotide complement of the sequence `seq`
"""
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

"""
`reverse(seq::NucleotideSequence)`

Reversed copy of the nucleotide sequence `seq`
"""
function reverse{T}(seq::NucleotideSequence{T})
    if isempty(seq)
        return seq
    end

    orphan!(seq)

    k = (2 * length(seq) + 63) % 64 + 1
    h = 64 - k

    data = zeros(Uint64, length(seq.data))
    j = length(data)
    i = 1
    @inbounds while true
        x = nucrev(seq.data[i])
        data[j] |= x >>> h
        if (j -= 1) == 0
            break
        end
        data[j] |= x << k;
        i += 1
    end

    return NucleotideSequence{T}(data, reverse(seq.ns), seq.part)
end

# Return the reverse complement of seq
# Return a reversed copy of seq
"""
`reverse_complement(seq::NucleotideSequence)`

Reversed complement of the nucleotide sequence `seq`
"""
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
    d, r = divrem64(i - 1)
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


# Find the index of the next 1 bit, starting at index i.
function nextone(b::BitVector, i)
    if i > length(b)
        return i
    end

    d, r = divrem64(i - 1)
    chunk = b.chunks[d + 1] >>> r
    if chunk != 0
        t = trailing_zeros(chunk)
        return i + t
    end

    i += 64 - r
    r = 0
    d += 1
    @inbounds for l in (d+1):length(b.chunks)
        if b.chunks[l] != 0
            return i + trailing_zeros(b.chunks[l])
        end
        i += 64
    end

    return min(length(b) + 1, i)
end


# Iterating through SequenceNIterator
# -----------------------------------

start(it::SequenceNIterator) = nextn(it, it.part.start)

function next(it::SequenceNIterator, i)
    d, r = divrem64(i - 1)
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

"Mismatch count between two kmers"
function nucmismatches(x::Uint64, y::Uint64)
    xyxor = x $ y
    return count_ones((xyxor & 0x5555555555555555) | ((xyxor & 0xAAAAAAAAAAAAAAAA) >>> 1))
end

"Mask of the first `k` bits of a Uint64"
function makemask(k::Integer)
    return 0xffffffffffffffff >> (64 - k)
end

"""
`mismatches(a::NucleotideSequence, b::NucleotideSequence, [nmatches=false])`

Return the number of mismatches between `a` and `b`.

If `a` and `b` are of differing lengths, only the first `min(length(a), length(b))`
nucleotides are compared.

# Arguments
* `a`: first sequence to compare
* `b`: second sequence to compare
* `nmatches`: if true, N matches anything, if false, N matches only itself (false)

# Returns
The number of mismatches
"""
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
    d0, r0 = divrem64(b.part.start - 1)
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

"Convert a Uint64 to a DNAKmer"
convert{K}(::Type{DNAKmer{K}}, x::Uint64) = box(DNAKmer{K}, unbox(Uint64, x))

"Convert a Uint64 to a RNAKmer"
convert{K}(::Type{RNAKmer{K}}, x::Uint64) = box(RNAKmer{K}, unbox(Uint64, x))

"Convert a DNAKmer to a Uint64"
convert{K}(::Type{Uint64}, x::DNAKmer{K})       = box(Uint64, unbox(DNAKmer{K}, x))

"Convert a RNAKmer to a Uint64"
convert{K}(::Type{Uint64}, x::RNAKmer{K})       = box(Uint64, unbox(RNAKmer{K}, x))


# Conversion to/from String

"Convert a String to a Kmer"
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

"Convert a String to a Kmer"
convert{T}(::Type{Kmer{T}}, seq::String) = convert(Kmer{T, length(seq)}, seq)

"Convert a Kmer to a String"
convert{T, K}(::Type{String}, seq::Kmer{T, K}) = convert(String, [convert(Char, x) for x in seq])


# Conversion to/from NucleotideSequence

"Convert a NucleotideSequence to a Kmer"
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

"Convert a NucleotideSequence to a Kmer"
convert{T}(::Type{Kmer}, seq::NucleotideSequence{T})    = convert(Kmer{T, length(seq)}, seq)

"Convert a NucleotideSequence to a Kmer"
convert{T}(::Type{Kmer{T}}, seq::NucleotideSequence{T}) = convert(Kmer{T, length(seq)}, seq)

"Convert a Kmer to a NucleotideSequence"
function convert{T, K}(::Type{NucleotideSequence{T}}, x::Kmer{T, K})
    ns = BitVector(K)
    fill!(ns, false)
    return NucleotideSequence{T}([convert(Uint64, x)], ns, 1:K)
end

"Convert a Kmer to a NucleotideSequence"
convert{T, K}(::Type{NucleotideSequence}, x::Kmer{T, K}) = convert(NucleotideSequence{T}, x)




# Constructors
# ------------

# From strings

"Construct a DNAKmer to a String"
dnakmer(seq::String) = convert(DNAKmer, seq)

"Construct a RNAKmer to a String"
rnakmer(seq::String) = convert(RNAKmer, seq)

"Construct a Kmer from a sequence of Nucleotides"
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

"Construct a Kmer from a DNASequence"
function kmer(seq::DNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(DNAKmer{length(seq)}, seq)
end

"Construct a Kmer from a RNASequence"
function kmer(seq::RNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(RNAKmer{length(seq)}, seq)
end

# call kmer with @inline macro would reduce the performance significantly?
# Would the compiler inline even without @inline?
"Construct a DNAKmer from a DNASequence"
function dnakmer(seq::DNASequence)
    @assert length(seq) <= 32 error("Cannot construct a K-mer longer than 32nt.")
    return convert(DNAKmer{length(seq)}, seq)
end

"Construct a RNAKmer from a RNASequence"
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

"""
`complement(kmer::Kmer)`

The Kmer complement of `kmer`
"""
function complement{T, K}(x::Kmer{T, K})
    return convert(Kmer{T, K},
        (~convert(Uint64, x)) & (0xffffffffffffffff >>> (2 * (32 - K))))
end

"""
`reverse(kmer::Kmer)`

Reversed copy of `kmer`
"""
reverse{T, K}(x::Kmer{T, K}) = convert(Kmer{T, K}, nucrev(convert(Uint64, x)) >>> (2 * (32 - K)))

"""
`reverse_complement(kmer::Kmer)`

Reversed complement of `kmer`
"""
reverse_complement{T, K}(x::Kmer{T, K}) = complement(reverse(x))

"""
`mismatches(x::Kmer, y::Kmer)`

Return the number of mismatches between `x` and `y`.

# Arguments
* `x`: first sequence to compare
* `y`: second sequence to compare

# Returns
The number of mismatches
"""
mismatches{T, K}(x::Kmer{T, K}, y::Kmer{T, K}) = nucmismatches(convert(Uint64, x), convert(Uint64, y))

"""
`canonical(x::Kmer)`

Canonical k-mer of `x`

A canonical k-mer is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting k-mers in data that is not strand specific,
and thus observing k-mer is equivalent to observing its reverse complement.
"""
function canonical{T, K}(x::Kmer{T, K})
    y = reverse_complement(x)
    return x < y ? x : y
end


"""
Iterate through k-mers neighboring on a de Bruijn graph.
"""
function neighbors{T, K}(x::Kmer{T, K})
    return KmerNeighborIterator{T, K}(x)
end


immutable KmerNeighborIterator{T, K}
    x::Kmer{T, K}
end


start(it::KmerNeighborIterator) = @compat UInt64(0)
done(it::KmerNeighborIterator, i) = i == 4
length(::KmerNeighborIterator) = 4
next{T, K}(it::KmerNeighborIterator{T, K}, i) =
    convert(Kmer{T, K}, (convert(Uint64, it.x) >>> 2) | (i << (2 * K - 2))), i + 1


# EachKmerIterator and EachKmerIteratorState
# ==========================================

# Iterate through every k-mer in a nucleotide sequence
immutable EachKmerIterator{T, K}
    seq::NucleotideSequence{T}
    step::Int
end


immutable EachKmerIteratorState{T, K}
    i::Int
    npos::Int
end

# Maybe this function should replace the default constructor.
# Is the (unsafe) default constructor used throughout our code?
"""
`eachkmer(seq::NucleotideSequence, k, [step=1])`

Construct a EachKmerIterator from a NucleotideSequence `seq` of size `k` and, optionally, a `step` value.

Differently from the default EachKmerIterator constructor, this function checks the validity of the arguments before construction.

# Arguments
* `seq`: A NucleotideSequence
* `k`: The size of each Kmer
* `step`: Number of steps

# Returns
A EachKmerIterator constructed with these parameters
"""
function each{T, K}(::Type{Kmer{T, K}}, seq::NucleotideSequence{T}, step::Integer=1)
    @assert K >= 0 "K must be ≥ 0 in EachKmer"
    @assert K <= 32 "K must be ≤ 32 in EachKmer"
    @assert step >= 1 "step must be ≥ 1"

    return EachKmerIterator{T, K}(seq, step)
end


function eachkmer{T}(seq::NucleotideSequence{T}, k::Integer, step::Integer=1)
    return each(Kmer{T, K}, seq, step)
end


@inline function getkmer{T, K}(::Type{Kmer{T, K}}, seq::NucleotideSequence{T}, i)
    d, r = divrem32(i - 1)
    x = seq.data[d + 1] >>> (2*r)
    if r + K > 32
        x |= seq.data[d + 2] << (2 * (32 - r))
    end
    mask = 0xffffffffffffffff >>> (64 - (2*K))
    return convert(Kmer{T, K}, x & mask)
end


@inline function nextkmerpos{T, K}(it::EachKmerIterator{T, K}, i, npos)
    while i + K - 1 >= npos
        while npos < i
            npos = nextone(it.seq.ns, i)
        end

        if i + K - 1 < npos
            break
        end

        i += it.step
        if i + K - 1 > it.seq.part.stop
            return EachKmerIteratorState{T, K}(i, npos)
        end
    end
    return EachKmerIteratorState{T, K}(i, npos)
end


function start{T, K}(it::EachKmerIterator{T, K})
    if K == 0
        return EachKmerIteratorState{T, K}(it.seq.part.stop + 2, 0)
    else
        npos = nextone(it.seq.ns, it.seq.part.start)
        return nextkmerpos(it, 1, npos)
    end
end


@inline function next{T, K}(it::EachKmerIterator{T, K},
                    state::EachKmerIteratorState{T, K})
    value = getkmer(Kmer{T, K}, it.seq, state.i)
    next_state = nextkmerpos(it, state.i + it.step, state.npos)
    return (state.i, value), next_state
end


@inline function done{T, K}(it::EachKmerIterator{T, K},
                            state::EachKmerIteratorState{T, K})
    return state.i + K - 1 > it.seq.part.stop
end


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

"""
`NucleotideCounts(seq::NucleotideSequence)`

Constructs a NucleotideCounts object from a NucleotideSequence `seq`.
"""
function NucleotideCounts{T}(seq::NucleotideSequence{T})
    dn, rn = divrem64(seq.part.start - 1)

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

"""
`NucleotideCounts(seq::Kmer)`

Constructs a NucleotideCounts object from a Kmer `seq`.
"""
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


# Kmer Composition
# ----------------

"""
Count occurances of short (<= 32) k-mers in a sequence.

# Arguments:
  * 'seq`: A NucleotideSequence
  * `step`: K-mers counted are separated by this many nucleotides (deafult: 1)
"""
immutable KmerCounts{T, K}
    data::Vector{Uint32}

    function KmerCounts(seq::NucleotideSequence{T}, step::Integer=1)
        data = zeros(Uint32, 4^K)
        @inbounds for (_, x) in each(Kmer{T, K}, seq, step)
            data[convert(Uint64, x) + 1] += 1
        end
        return new(data)
    end
end

typealias DNAKmerCounts{K} KmerCounts{DNANucleotide, K}
typealias RNAKmerCounts{K} KmerCounts{DNANucleotide, K}


function getindex{T, K}(counts::KmerCounts{T, K}, x::Kmer{T, K})
    @inbounds c = counts.data[convert(Uint64, x) + 1]
    return c
end


function show{T, K}(io::IO, counts::KmerCounts{T, K})
    println(io, (T == DNANucleotide ? "DNA" : "RNA"), "KmerCounts{", K, "}:")
    for x in (@compat UInt64(1)):(@compat UInt64(4^K))
        s = convert(String, convert(Kmer{T, K}, x - 1))
        println(io, "  ", s, " => ", counts.data[x])
    end
end


