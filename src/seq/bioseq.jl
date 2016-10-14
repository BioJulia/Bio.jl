# BioSequence
# ===========
#
# A general purpose biological sequence representation.
#
# It is cheap to create a subsequence from a sequence because sequences can
# share underlying data: creating a subsequence is just allocating a new
# sequence object defining a part of the original sequence without copy.
# Destructive operations create a copy of the underlying data if and only if it
# is shared between (sub)sequences. This is often called as copy-on-write
# strategy in computer science and should be transparent to the user.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# About Internals
# ---------------
#
# The `data` field of a `BioSequence{A}` object contains binary representation
# of a biological character sequence. Each character is encoded with an encoder
# corresponding to the alphabet `A` and compactly packed into `data`. To extract
# a character from a sequence, you should decode this binary sequence with a
# decoder that is a pair of the encoder. The length of encoded binary bits is
# fixed, and hence a character at arbitrary position can be extracted in a
# constant time. To know the exact location of a character at a position, you
# can use the `bitsof(seq, i)` function, which returns a pair of element's index
# containing binary bits and bits' offset. As a whole, character extraction
# `seq[i]` can be written as:
#
#     j = bitindex(seq, i)
#     decode(A, (seq.data[index(j)] >> offset(j)) & mask(A))
#
#  index :        index(j) - 1       index(j)       index(j) + 1
#   data : ....|xxxxx...........|xxXxxxxxxxxxxxxx|............xxxx|....
# offset :                          |<-offset(j)-|
#  width :      |<---- 64 ---->| |<---- 64 ---->| |<---- 64 ---->|
#
#  * '.' : unused (4 bits/char)
#  * 'x' : used
#  * 'X' : used and pointed by index `i`

"""
Biological sequence data structure indexed by an alphabet type `A`.
"""
type BioSequence{A<:Alphabet} <: Sequence
    data::Vector{UInt64}  # encoded character sequence data
    part::UnitRange{Int}  # interval within `data` defining the (sub)sequence
    shared::Bool          # true if and only if `data` is shared between sequences

    function BioSequence(data::Vector{UInt64}, part::UnitRange{Int}, shared::Bool)
        return new(data, part, shared)
    end
end

typealias DNASequence       BioSequence{DNAAlphabet{4}}
typealias RNASequence       BioSequence{RNAAlphabet{4}}
typealias AminoAcidSequence BioSequence{AminoAcidAlphabet}
typealias CharSequence      BioSequence{CharAlphabet}


# String Decorators
# -----------------

remove_newlines(s) = replace(s, r"\r|\n", "")

macro dna_str(seq, flag)
    if flag == "s"
        return DNASequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            DNASequence($(remove_newlines(seq)))
        end
    end
    error("Invalid DNA flag: '$(flag)'")
end

macro dna_str(seq)
    return DNASequence(remove_newlines(seq))
end

macro rna_str(seq, flag)
    if flag == "s"
        return RNASequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            RNASequence($(remove_newlines(seq)))
        end
    end
    error("Invalid RNA flag: '$(flag)'")
end

macro rna_str(seq)
    return RNASequence(remove_newlines(seq))
end

macro aa_str(seq, flag)
    if flag == "s"
        return AminoAcidSequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            AminoAcidSequence($(remove_newlines(seq)))
        end
    end
    error("Invalid Amino Acid flag: '$(flag)'")
end

macro aa_str(seq)
    return AminoAcidSequence(remove_newlines(seq))
end

macro char_str(seq, flag)
    if flag == "s"
        return CharSequence(remove_newlines(seq))
    elseif flag == "d"
        return quote
            CharSequence($(remove_newlines(seq)))
        end
    end
    error("Invalid Char flag: '$(flag)'")
end

macro char_str(seq)
    return CharSequence(remove_newlines(seq))
end


# Constructors
# ------------

function seq_data_len{A}(::Type{A}, len::Integer)
    return cld(len, div(64, bitsof(A)))
end

function (::Type{BioSequence{A}}){A<:Alphabet}(len::Integer)
    return BioSequence{A}(Vector{UInt64}(seq_data_len(A, len)), 1:len, false)
end

BioSequence(::Type{DNANucleotide}) = DNASequence()
BioSequence(::Type{RNANucleotide}) = RNASequence()
BioSequence(::Type{AminoAcid}) = AminoAcidSequence()
BioSequence(::Type{Char}) = CharSequence()

function BioSequence()
    return BioSequence{VoidAlphabet}(Vector{UInt64}(0), 0:-1, false)
end

function (::Type{BioSequence{A}}){A<:Alphabet}(
        src::Union{AbstractString,AbstractVector},
        startpos::Integer=1,
        stoppos::Integer=length(src))
    len = stoppos - startpos + 1
    seq = BioSequence{A}(len)
    return encode_copy!(seq, 1, src, startpos, len)
end

# create a subsequence
function BioSequence{A,T<:Integer}(other::BioSequence{A}, part::UnitRange{T})
    checkbounds(other, part)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = BioSequence{A}(other.data, start:stop, true)
    other.shared = true
    return subseq
end

function (::Type{BioSequence{A}}){A}(other::BioSequence{A}, part::UnitRange)
    return BioSequence(other, part)
end

# concatenate chunks
function (::Type{BioSequence{A}}){A}(chunks::BioSequence{A}...)
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seq = BioSequence{A}(len)
    offset = 1
    for chunk in chunks
        copy!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

function Base.repeat{A}(chunk::BioSequence{A}, n::Integer)
    seq = BioSequence{A}(length(chunk) * n)
    offset = 1
    for _ in 1:n
        copy!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

# operators for concat and repeat
Base.:*{A}(chunk::BioSequence{A}, chunks::BioSequence{A}...) =
    BioSequence{A}(chunk, chunks...)
Base.:^(chunk::BioSequence, n::Integer) = repeat(chunk, n)

# conversion between different alphabet size
for A in [DNAAlphabet, RNAAlphabet]
    # 4 bits => 2 bits
    @eval function Base.convert(::Type{BioSequence{$(A{2})}}, seq::BioSequence{$(A{4})})
        # TODO: make it faster with bit-parallel algorithm
        newseq = BioSequence{$(A{2})}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end

    # 2 bits => 4 bits
    @eval function Base.convert(::Type{BioSequence{$(A{4})}}, seq::BioSequence{$(A{2})})
        newseq = BioSequence{$(A{4})}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end
end

# conversion between DNA and RNA
for (A1, A2) in [(DNAAlphabet, RNAAlphabet), (RNAAlphabet, DNAAlphabet)], n in (2, 4)
    # NOTE: assumes that binary representation is identical between DNA and RNA
    @eval function Base.convert(::Type{BioSequence{$(A1{n})}},
                                seq::BioSequence{$(A2{n})})
        newseq = BioSequence{$(A1{n})}(seq.data, seq.part, true)
        seq.shared = true
        return newseq
    end
end

# from a vector
function Base.convert{A<:DNAAlphabet}(::Type{BioSequence{A}},
                                      seq::AbstractVector{DNANucleotide})
    return BioSequence{A}(seq, 1, endof(seq))
end
function Base.convert{A<:RNAAlphabet}(::Type{BioSequence{A}},
                                      seq::AbstractVector{RNANucleotide})
    return BioSequence{A}(seq, 1, endof(seq))
end
function Base.convert(::Type{AminoAcidSequence}, seq::AbstractVector{AminoAcid})
    return AminoAcidSequence(seq, 1, endof(seq))
end

# to a vector
Base.convert(::Type{Vector}, seq::BioSequence) = collect(seq)
function Base.convert{A<:DNAAlphabet}(::Type{Vector{DNANucleotide}},
                                      seq::BioSequence{A})
    return collect(seq)
end
function Base.convert{A<:RNAAlphabet}(::Type{Vector{RNANucleotide}},
                                      seq::BioSequence{A})
    return collect(seq)
end
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = collect(seq)

# from/to a string
function Base.convert{S<:AbstractString}(::Type{S}, seq::BioSequence)
    return convert(S, [Char(x) for x in seq])
end
Base.convert{S<:AbstractString,A}(::Type{BioSequence{A}}, seq::S) = BioSequence{A}(seq)

# Summaries
# ---------

Base.summary{A<:DNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "DNA Sequence")
Base.summary{A<:RNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "RNA Sequence")
Base.summary(seq::AminoAcidSequence) = string(length(seq), "aa ", "Amino Acid Sequence")
Base.summary(seq::CharSequence) = string(length(seq), "char ", "Char Sequence")


# Basic Operators
# ---------------

alphabet{A}(::Type{BioSequence{A}}) = alphabet(A)

Base.length(seq::BioSequence) = length(seq.part)
Base.eltype{A}(::Type{BioSequence{A}}) = eltype(A)

function Base.count(f::Function, seq::BioSequence)
    n = 0
    for x in seq
        if f(x)
            n += 1
        end
    end
    return n
end

function Base.map(f::Function, seq::BioSequence)
    return map!(f, copy(seq))
end

function Base.filter(f::Function, seq::BioSequence)
    return filter!(f, copy(seq))
end

function Base.checkbounds(seq::BioSequence, range::UnitRange)
    if 1 ≤ range.start && range.stop ≤ endof(seq)
        return true
    end
    throw(BoundsError(seq, range))
end

function Base.checkbounds(seq::BioSequence, locs::AbstractVector{Bool})
    if length(seq) == length(locs)
        return true
    end
    throw(BoundsError(seq, locs))
end

function Base.checkbounds(seq::BioSequence, locs::AbstractVector)
    for i in locs
        checkbounds(seq, i)
    end
    return true
end

function checkdimension(from::Integer, to::Integer)
    if from == to
        return true
    end
    throw(DimensionMismatch(string(
        "attempt to assign ",
        from, " elements to ",
        to,   " elements")))
end

function checkdimension(seq::BioSequence, locs::AbstractVector)
    return checkdimension(length(seq), length(locs))
end

function checkdimension(seq::BioSequence, locs::AbstractVector{Bool})
    return checkdimension(length(seq), sum(locs))
end

# creates a bit mask for given number of bits `n`
mask(n::Integer) = (UInt64(1) << n) - 1
mask{A<:Alphabet}(::Type{A}) = mask(bitsof(A))

# assumes `i` is positive and `bitsof(A)` is a power of 2
@inline function bitindex{A}(seq::BioSequence{A}, i::Integer)
    return BitIndex((i + first(seq.part) - 2) << trailing_zeros(bitsof(A)))
end

@inline function inbounds_getindex{A}(seq::BioSequence{A}, i::Integer)
    j = bitindex(seq, i)
    @inbounds return decode(A, (seq.data[index(j)] >> offset(j)) & mask(A))
end

Base.getindex(seq::BioSequence, part::UnitRange) = BioSequence(seq, part)
Base.view(seq::BioSequence, part::UnitRange) = BioSequence(seq, part)

function Base.setindex!{A,T<:Integer}(seq::BioSequence{A},
                                      other::BioSequence{A},
                                      locs::AbstractVector{T})
    checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    for (i, x) in zip(locs, other)
        unsafe_setindex!(seq, x, i)
    end
    return seq
end

function Base.setindex!{A,T<:Integer}(seq::BioSequence{A},
                                      other::BioSequence{A},
                                      locs::UnitRange{T})
    checkbounds(seq, locs)
    checkdimension(other, locs)
    return copy!(seq, locs.start, other, 1)
end

function Base.setindex!{A}(seq::BioSequence{A},
                           other::BioSequence{A},
                           locs::AbstractVector{Bool})
    checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    i = j = 0
    while (i = findnext(locs, i + 1)) > 0
        unsafe_setindex!(seq, other[j+=1], i)
    end
    return seq
end

function Base.setindex!{A}(seq::BioSequence{A}, other::BioSequence{A}, ::Colon)
    return setindex!(seq, other, 1:endof(seq))
end

function Base.setindex!(seq::BioSequence, x, i::Integer)
    checkbounds(seq, i)
    orphan!(seq)
    return unsafe_setindex!(seq, x, i)
end

function Base.setindex!{A,T<:Integer}(seq::BioSequence{A}, x, locs::AbstractVector{T})
    checkbounds(seq, locs)
    bin = enc64(seq, x)
    orphan!(seq)
    for i in locs
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

function Base.setindex!{A}(seq::BioSequence{A}, x, locs::AbstractVector{Bool})
    checkbounds(seq, locs)
    bin = enc64(seq, x)
    orphan!(seq)
    i = j = 0
    while (i = findnext(locs, i + 1)) > 0
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

Base.setindex!{A}(seq::BioSequence{A}, x, ::Colon) = setindex!(seq, x, 1:endof(seq))

# this is "unsafe" because of no bounds check and no orphan! call
@inline function unsafe_setindex!{A}(seq::BioSequence{A}, x, i::Integer)
    bin = enc64(seq, x)
    return encoded_setindex!(seq, bin, i)
end

@inline function encoded_setindex!{A}(seq::BioSequence{A}, bin::UInt64, i::Integer)
    j, r = bitindex(seq, i)
    @inbounds seq.data[j] = (bin << r) | (seq.data[j] & ~(mask(A) << r))
    return seq
end

function Base.resize!{A}(seq::BioSequence{A}, size::Integer)
    if size < 0
        throw(ArgumentError("size must be non-negative"))
    end
    orphan!(seq, size)
    resize!(seq.data, seq_data_len(A, size + seq.part.start - 1))
    seq.part = seq.part.start:seq.part.start+size-1
    return seq
end

Base.empty!(seq::BioSequence) = resize!(seq, 0)

function Base.push!{A}(seq::BioSequence{A}, x)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    encoded_setindex!(seq, bin, endof(seq))
    return seq
end

function Base.unshift!{A}(seq::BioSequence{A}, x)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    copy!(seq, 2, seq, 1, length(seq) - 1)
    encoded_setindex!(seq, bin, 1)
    return seq
end

function Base.pop!(seq::BioSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    x = seq[end]
    deleteat!(seq, endof(seq))
    return x
end

function Base.shift!(seq::BioSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    x = seq[1]
    deleteat!(seq, 1)
    return x
end

function Base.insert!{A}(seq::BioSequence{A}, i::Integer, x)
    checkbounds(seq, i)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    copy!(seq, i + 1, seq, i, endof(seq) - i)
    encoded_setindex!(seq, bin, i)
    return seq
end

function Base.deleteat!{A,T<:Integer}(seq::BioSequence{A}, range::UnitRange{T})
    checkbounds(seq, range)
    copy!(seq, range.start, seq, range.stop + 1, length(seq) - range.stop)
    resize!(seq, length(seq) - length(range))
    return seq
end

function Base.deleteat!(seq::BioSequence, i::Integer)
    checkbounds(seq, i)
    copy!(seq, i, seq, i + 1, length(seq) - i)
    resize!(seq, length(seq) - 1)
    return seq
end

function Base.append!{A}(seq::BioSequence{A}, other::BioSequence{A})
    resize!(seq, length(seq) + length(other))
    copy!(seq, endof(seq) - length(other) + 1, other, 1)
    return seq
end

function Base.copy!{A}(seq::BioSequence{A}, doff::Integer,
                       src::Vector{UInt8},  soff::Integer, len::Integer)
    datalen = seq_data_len(A, len)
    return seq
end

function Base.copy!{A}(dst::BioSequence{A}, src::BioSequence{A})
    return copy!(dst, 1, src, 1)
end

function Base.copy!{A}(dst::BioSequence{A}, doff::Integer,
                       src::BioSequence{A}, soff::Integer)
    return copy!(dst, doff, src, soff, length(src) - soff + 1)
end

function Base.copy!{A}(dst::BioSequence{A}, doff::Integer,
                       src::BioSequence{A}, soff::Integer, len::Integer)
    checkbounds(dst, doff:doff+len-1)
    checkbounds(src, soff:soff+len-1)

    if dst.shared || (dst === src && doff > soff)
        orphan!(dst, length(dst), true)
    end

    id = bitindex(dst, doff)
    is = bitindex(src, soff)
    rest = len * bitsof(A)

    while rest > 0
        # move `k` bits from `src` to `dst`
        x = dst.data[index(id)]
        y = src.data[index(is)]
        if offset(id) < offset(is)
            y >>= offset(is) - offset(id)
            k = min(64 - offset(is), rest)
        else
            y <<= offset(id) - offset(is)
            k = min(64 - offset(id), rest)
        end
        m = mask(k) << offset(id)
        dst.data[index(id)] = y & m | x & ~m

        id += k
        is += k
        rest -= k
    end

    return dst
end

function Base.map!(f::Function, seq::BioSequence)
    orphan!(seq)
    for i in 1:endof(seq)
        unsafe_setindex!(seq, f(inbounds_getindex(seq, i)), i)
    end
    return seq
end

function Base.filter!{A}(f::Function, seq::BioSequence{A})
    orphan!(seq)

    len = 0
    next = bitindex(seq, 1)
    j = index(next)
    datum::UInt64 = 0
    for i in 1:endof(seq)
        x = inbounds_getindex(seq, i)
        if f(x)
            datum |= enc64(seq, x) << offset(next)
            len += 1
            next += bitsof(A)
            if index(next) != j
                seq.data[j] = datum
                datum = 0
                j = index(next)
            end
        end
    end
    if offset(next) > 0
        seq.data[j] = datum
    end
    resize!(seq, len)

    return seq
end

function Base.similar{A}(seq::BioSequence{A}, len::Integer=length(seq))
    return BioSequence{A}(len)
end

# actually, users don't need to create a copy of a sequence.
function Base.copy{A}(seq::BioSequence{A})
    newseq = BioSequence{A}(seq, 1:endof(seq))
    orphan!(newseq, length(seq), true)  # force orphan!
    @assert newseq.data !== seq.data
    return newseq
end

# Replace a BioSequence's data with a copy, copying only what's needed.
# The user should never need to call this, as it has no outward effect on the
# sequence.
function orphan!{A}(seq::BioSequence{A}, size::Integer=length(seq), force::Bool=false)
    if !seq.shared && !force
        return seq
    end

    j, r = bitindex(seq, 1)
    data = Vector{UInt64}(seq_data_len(A, size))

    if !isempty(seq) && !isempty(data)
        x = seq.data[j] >> r
        m = index(bitindex(seq, endof(seq))) - j + 1
        l = min(endof(data), m)
        @inbounds @simd for i in 1:l-1
            y = seq.data[j+i]
            data[i] = x | y << (64 - r)
            x = y >> r
        end
        if m <= l
            data[l] = x
        else
            y = seq.data[j+l]
            data[l] = x | y << (64 - r)
        end
    end

    seq.data = data
    seq.part = 1:length(seq)
    seq.shared = false
    return seq
end


# Transformations
# ---------------

function Base.reverse!(seq::BioSequence)
    orphan!(seq)
    for i in 1:div(endof(seq), 2)
        x = inbounds_getindex(seq, i)
        unsafe_setindex!(seq, inbounds_getindex(seq, endof(seq) - i + 1), i)
        unsafe_setindex!(seq, x, endof(seq) - i + 1)
    end
    return seq
end

Base.reverse(seq::BioSequence) = reverse!(copy(seq))

@generated function Base.reverse{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    n = bitsof(A)
    if n == 2
        nucrev = :nucrev2
    elseif n == 4
        nucrev = :nucrev4
    else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        data = Vector{UInt64}(seq_data_len(A, length(seq)))
        i = 1
        next = bitindex(seq, endof(seq))
        stop = bitindex(seq, 0)
        r = rem(offset(next) + $n, 64)
        if r == 0
            @inbounds while next - stop > 0
                x = seq.data[index(next)]
                data[i] = $nucrev(x)
                i += 1
                next -= 64
            end
        else
            @inbounds while next - stop > 64
                j = index(next)
                x = (seq.data[j] << (64 - r)) | (seq.data[j-1] >> r)
                data[i] = $nucrev(x)
                i += 1
                next -= 64
            end
            if next - stop > 0
                data[i] = $nucrev(seq.data[index(next)] << (64 - r))
            end
        end
        return BioSequence{A}(data, 1:length(seq), false)
    end
end

@inline function nucrev2(x::UInt64)
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>  2
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
     return x
end

@inline function nucrev4(x::UInt64)
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
     return x
end

"""
    complement!(seq)

Make a complement sequence of `seq` in place.
"""
function complement!{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(seq::BioSequence{A})
    orphan!(seq)
    next = bitindex(seq, 1)
    stop = bitindex(seq, endof(seq) + 1)
    @inbounds while next < stop
        seq.data[index(next)] = ~seq.data[index(next)]
        next += 64
    end
    return seq
end

function complement!{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(seq::BioSequence{A})
    orphan!(seq)
    next = bitindex(seq, 1)
    stop = bitindex(seq, endof(seq) + 1)
    @inbounds while next < stop
        x = seq.data[index(next)]
        seq.data[index(next)] = (
            ((x & 0x1111111111111111) << 3) | ((x & 0x8888888888888888) >> 3) |
            ((x & 0x2222222222222222) << 1) | ((x & 0x4444444444444444) >> 1))
        next += 64
    end
    return seq
end

"""
    complement(seq)

Make a complement sequence of `seq`.
"""
function complement{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(copy(seq))
end

"""
    reverse_complement!(seq)

Make a reversed complement sequence of `seq` in place.

Ambiguous nucleotides are left as-is.
"""
function reverse_complement!{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(reverse!(seq))
end

"""
    reverse_complement(seq)

Make a reversed complement sequence of `seq`.

Ambiguous nucleotides are left as-is.
"""
function reverse_complement{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(reverse(seq))
end


# Ambiguous nucleotides iterator
# ------------------------------

immutable AmbiguousNucleotideIterator{A<:Union{DNAAlphabet,RNAAlphabet}}
    seq::BioSequence{A}
end

ambiguous_positions(seq::BioSequence) = AmbiguousNucleotideIterator(seq)

Base.start(it::AmbiguousNucleotideIterator) = find_next_ambiguous(it.seq, 1)
Base.done(it::AmbiguousNucleotideIterator, nextpos) = nextpos == 0
function Base.next(it::AmbiguousNucleotideIterator, nextpos)
    return nextpos, find_next_ambiguous(it.seq, nextpos + 1)
end

Base.iteratorsize(::AmbiguousNucleotideIterator) = Base.SizeUnknown()

function find_next_ambiguous{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(
        seq::BioSequence{A}, i::Integer)
    # no ambiguity
    return 0
end

function find_next_ambiguous{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(
        seq::BioSequence{A}, from::Integer)
    for i in max(from, 1):endof(seq)
        nt = inbounds_getindex(seq, i)
        if isambiguous(nt)
            return i
        end
    end
    # no ambiguity
    return 0
end


# Mismatch counting
# -----------------

"""
    mismatches(seq1::BioSequence, seq2::BioSequence[, compatible=false])

Return the number of mismatches between `seq1` and `seq2`.

If `seq1` and `seq2` are of differing lengths, only the first `min(length(seq1),
length(seq2))` nucleotides are compared.  When `compatible` is `true`, sequence
symbols are comapred using `iscompatible`; otherwise using `==`.
"""
function mismatches{A<:Alphabet}(
        seq1::BioSequence{A},
        seq2::BioSequence{A},
        compatible::Bool=false)
    if ((bitsof(A) == 2 || bitsof(A) == 4) && !compatible) ||
        A == DNAAlphabet{2} ||
        A == RNAAlphabet{2}
        return bitparallel_mismatches(seq1, seq2)
    end

    mis = 0
    if compatible
        for (x, y) in zip(seq1, seq2)
            if !iscompatible(x, y)
                mis += 1
            end
        end
    else
        for (x, y) in zip(seq1, seq2)
            if x != y
                mis += 1
            end
        end
    end
    return mis
end

@generated function bitparallel_mismatches{A}(a::BioSequence{A}, b::BioSequence{A})
    n = bitsof(A)
    if n == 2
        bitpar_mismatches = :bitpar_mismatches2
    elseif n == 4
        bitpar_mismatches = :bitpar_mismatches4
    else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        if length(a) > length(b)
            return bitparallel_mismatches(b, a)
        end
        @assert length(a) ≤ length(b)

        nexta = bitindex(a, 1)
        nextb = bitindex(b, 1)
        stopa = bitindex(a, endof(a) + 1)
        mismatches = 0

        # align reading position of `a.data` so that `offset(nexta) == 0`
        if nexta < stopa && offset(nexta) != 0
            x = a.data[index(nexta)] >> offset(nexta)
            y = b.data[index(nextb)] >> offset(nextb)
            if offset(nextb) > offset(nexta)
                y |= b.data[index(nextb)+1] << (64 - offset(nextb))
            end
            k = 64 - offset(nexta)
            m = mask(k)
            mismatches += $bitpar_mismatches(x & m, y & m)
            nexta += k
            nextb += k
        end
        @assert offset(nexta) == 0

        if offset(nextb) == 0  # data are aligned with each other
            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                mismatches += $bitpar_mismatches(x, y)
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                m = mask(stopa - nexta)
                mismatches += $bitpar_mismatches(x & m, y & m)
            end
        elseif nexta < stopa
            y = b.data[index(nextb)]
            nextb += 64

            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                z = b.data[index(nextb)]
                y = y >> offset(nextb) | z << (64 - offset(nextb))
                mismatches += $bitpar_mismatches(x, y)
                y = z
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                x = a.data[index(nexta)]
                y = y >> offset(nextb)
                if 64 - offset(nextb) < stopa - nexta
                    y |= a.data[index(nextb)] << (64 - offset(nextb))
                end
                m = mask(stopa - nexta)
                mismatches += $bitpar_mismatches(x & m, y & m)
            end
        end

        return mismatches
    end
end

# bit-parallel mismatch count algorithm for 2 and 4-bit encoding
@inline function bitpar_mismatches2(x::UInt64, y::UInt64)
    xyxor = x $ y
    mismatches = UInt64(0)
    mismatches |=  xyxor & 0x5555555555555555
    mismatches |= (xyxor & 0xAAAAAAAAAAAAAAAA) >> 1
    return count_ones(mismatches)
end

@inline function bitpar_mismatches4(x::UInt64, y::UInt64)
    xyxor = x $ y
    mismatches = UInt64(0)
    mismatches |=  xyxor & 0x1111111111111111
    mismatches |= (xyxor & 0x2222222222222222) >> 1
    mismatches |= (xyxor & 0x4444444444444444) >> 2
    mismatches |= (xyxor & 0x8888888888888888) >> 3
    return count_ones(mismatches)
end


# Shuffle
# -------

function Base.shuffle(seq::BioSequence)
    return shuffle!(copy(seq))
end

function Base.shuffle!(seq::BioSequence)
    orphan!(seq)
    # Fisher-Yates shuffle
    for i in 1:endof(seq)-1
        j = rand(i:endof(seq))
        seq[i], seq[j] = seq[j], seq[i]
    end
    return seq
end


# Encoding
# --------

function encode_copy!{A}(dst::BioSequence{A},
                         src::Union{AbstractVector,AbstractString})
    return encode_copy!(dst, 1, src, 1)
end

function encode_copy!{A}(dst::BioSequence{A},
                         doff::Integer,
                         src::Union{AbstractVector,AbstractString},
                         soff::Integer)
    return encode_copy!(dst, doff, src, soff, length(src) - soff + 1)
end

function encode_copy!{A}(dst::BioSequence{A},
                         doff::Integer,
                         src::Union{AbstractVector,AbstractString},
                         soff::Integer,
                         len::Integer)
    if soff != 1 && !isascii(src)
        throw(ArgumentError("source offset ≠ 1 is not supported for non-ASCII string"))
    end

    checkbounds(dst, doff:doff+len-1)
    if length(src) < soff + len - 1
        throw(ArgumentError("source string does not contain $len elements from $soff"))
    end

    orphan!(dst)
    next = bitindex(dst, doff)
    stop = bitindex(dst, doff + len)
    i = soff
    while next < stop
        x = UInt64(0)
        j = index(next)
        while index(next) == j && next < stop
            char, i = Base.next(src, i)
            x |= enc64(dst, convert(Char, char)) << offset(next)
            next += bitsof(A)
        end
        dst.data[j] = x
    end
    return dst
end

function encode_copy!{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(
        dst::BioSequence{A}, doff::Integer,
        src::AbstractVector{UInt8}, soff::Integer, len::Integer)
    checkbounds(dst, doff:doff+len-1)
    if length(src) < soff + len - 1
        throw(ArgumentError("source string does not contain $len elements from $soff"))
    end

    orphan!(dst)
    charmap = A <: DNAAlphabet ? char_to_dna : char_to_rna
    i = soff
    next = bitindex(dst, doff)
    stop = bitindex(dst, doff + len)

    # head
    if offset(next) != 0
        for d in 0:div(64 - offset(next), 4)-1
            dst[doff+d] = charmap[src[i+d]+1]
        end
        i += div(64 - offset(next), 4)
        next += 64 - offset(next)
    end

    # body
    D = 16
    while next < (stop - offset(stop))
        x::UInt64 = 0
        check = 0x00
        @inbounds for d in 0:D-1
            y = reinterpret(UInt8, charmap[src[i+d]+1])
            x |= UInt64(y) << 4d
            check |= y
        end
        if check & 0x80 != 0
            # invalid byte(s) is detected
            for d in 0:D-1
                if !isvalid(charmap[src[i+d]+1])
                    error("cannot encode $(src[i+d])")
                end
            end
        end
        dst.data[index(next)] = x
        i += D
        next += 64
    end

    # tail
    for d in 0:div(stop - next, 4)-1
        dst[doff+i-soff+d] = charmap[src[i+d]+1]
    end

    return dst
end

function enc64{A}(::BioSequence{A}, x)
    return UInt64(encode(A, convert(eltype(A), x)))
end


# Sequences to Matrix
# -------------------

"""
    seqmatrix{A<:Alphabet}(vseq::Vector{BioSequence{A}}, major::Symbol)

Construct a matrix of nucleotides or amino acids from a vector of `BioSequence`s.

If parameter `major` is set to `:site`, the matrix is created such that one
nucleotide from each sequence is placed in each column i.e. the matrix is laid
out in site-major order.
This means that iteration over one position of many sequences is efficient,
as julia arrays are laid out in column major order.

If the parameter `major` is set to `:seq`, the matrix is created such that each
sequence is placed in one column i.e. the matrix is laid out in sequence-major
order.
This means that iteration across each sequence in turn is efficient,
as julia arrays are laid out in column major order.

# Examples
```julia
julia> seqs = [dna"AAA", dna"TTT", dna"CCC", dna"GGG"]
4-element Array{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},1}:
 3nt DNA Sequence:
AAA
 3nt DNA Sequence:
TTT
 3nt DNA Sequence:
CCC
 3nt DNA Sequence:
GGG

julia> seqmatrix(seqs, :site)
4x3 Array{Bio.Seq.DNANucleotide,2}:
 DNA_A  DNA_A  DNA_A
 DNA_T  DNA_T  DNA_T
 DNA_C  DNA_C  DNA_C
 DNA_G  DNA_G  DNA_G

 julia> seqmatrix(seqs, :seq)
 3x4 Array{Bio.Seq.DNANucleotide,2}:
  DNA_A  DNA_T  DNA_C  DNA_G
  DNA_A  DNA_T  DNA_C  DNA_G
  DNA_A  DNA_T  DNA_C  DNA_G
```
"""
function seqmatrix{A<:Alphabet}(vseq::AbstractVector{BioSequence{A}}, major::Symbol)
    nseqs = length(vseq)
    @assert nseqs > 0 throw(ArgumentError("Vector of BioSequence{$A} is empty."))
    nsites = length(vseq[1])
    @inbounds for i in 2:nseqs
        length(vseq[i]) == nsites || throw(ArgumentError("Sequences in vseq must be of same length"))
    end
    if major == :site
        mat = Matrix{eltype(A)}(nseqs, nsites)
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[seq, site] = vseq[seq][site]
        end
        return mat
    elseif major == :seq
        mat = Matrix{eltype(A)}(nsites, nseqs)
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[site, seq] = vseq[seq][site]
        end
        return mat
    else
        throw(ArgumentError("major must be :site or :seq"))
    end
end

"""
    seqmatrix{A<:Alphabet,T}(::Type{T}, vseq::Vector{BioSequence{A}}, major::Symbol)

Construct a matrix of `T` from a vector of `BioSequence`s.

If parameter `major` is set to `:site`, the matrix is created such that one
nucleotide from each sequence is placed in each column i.e. the matrix is laid
out in site-major order.
This means that iteration over one position of many sequences is efficient,
as julia arrays are laid out in column major order.

If the parameter `major` is set to `:seq`, the matrix is created such that each
sequence is placed in one column i.e. the matrix is laid out in sequence-major
order.
This means that iteration across each sequence in turn is efficient,
as julia arrays are laid out in column major order.

# Examples
```julia
julia> seqs = [dna"AAA", dna"TTT", dna"CCC", dna"GGG"]
4-element Array{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},1}:
 3nt DNA Sequence:
AAA
 3nt DNA Sequence:
TTT
 3nt DNA Sequence:
CCC
 3nt DNA Sequence:
GGG

julia> seqmatrix(seqs, :site, UInt8)
4×3 Array{UInt8,2}:
 0x01  0x01  0x01
 0x08  0x08  0x08
 0x02  0x02  0x02
 0x04  0x04  0x04

julia> seqmatrix(seqs, :seq, UInt8)
3×4 Array{UInt8,2}:
 0x01  0x08  0x02  0x04
 0x01  0x08  0x02  0x04
 0x01  0x08  0x02  0x04
```
"""
function seqmatrix{T,A<:Alphabet}(::Type{T}, vseq::AbstractVector{BioSequence{A}}, major::Symbol)
    nseqs = length(vseq)
    @assert nseqs > 0 throw(ArgumentError("Vector of BioSequence{$A} is empty."))
    nsites = length(vseq[1])
    @inbounds for i in 2:nseqs
        length(vseq[i]) == nsites || throw(ArgumentError("Sequences in vseq must be of same length."))
    end
    if major == :site
        mat = Matrix{T}(nseqs, nsites)
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[seq, site] = convert(T, vseq[seq][site])
        end
        return mat
    elseif major == :seq
        mat = Matrix{T}(nsites, nseqs)
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[site, seq] = convert(T, vseq[seq][site])
        end
        return mat
    else
        throw(ArgumentError("major must be :site or :seq"))
    end
end

# Consensus
# ---------

"""
    majorityvote{A<:NucleotideAlphabet}(seqs::AbstractVector{BioSequence{A}})

Construct a sequence that is a consensus of a vector of sequences.

The consensus is established by a simple majority vote rule, where amiguous
nucleotides cast an equal vote for each of their possible states.
For each site a winner(s) out of A, T(U), C, or G is determined, in the cases
of ties the ambiguity symbol that unifies all the winners is returned.
E.g if A and T tie, then W is inserted in the consensus. If all A, T, C, and G
tie at a site, then N is inserted in the consensus. Note this means that if a
nucletide e.g. 'C' and a gap '-' draw, the nucleotide will always win over the
gap, even though they tied.

# Examples
```julia
julia> seqs = [dna"CTCGATCGATCC", dna"CTCGAAAAATCA", dna"ATCGAAAAATCG", dna"ATCGGGGGATCG"]

4-element Array{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},1}:
 CTCGATCGATCC
 CTCGAAAAATCA
 ATCGAAAAATCG
 ATCGGGGGATCG

julia> majorityvote(seqs)
12nt DNA Sequence:
MTCGAAARATCG
```
"""
function majorityvote{A<:NucleotideAlphabet}(seqs::AbstractVector{BioSequence{A}})
    mat = seqmatrix(UInt8, seqs, :site)
    nsites = size(mat, 2)
    nseqs = size(mat, 1)
    result = BioSequence{A}(nsites)
    votes = Array{Int}(16)
    @inbounds for site in 1:nsites
        fill!(votes, 0)
        for seq in 1:nseqs
            nuc = mat[seq, site]
            votes[1] += nuc == 0x00
            votes[2] += (nuc & 0x01) != 0x00
            votes[3] += (nuc & 0x02) != 0x00
            votes[5] += (nuc & 0x04) != 0x00
            votes[9] += (nuc & 0x08) != 0x00
        end
        m = maximum(votes)
        merged = 0x00
        for i in 0x01:0x10
            merged |= ifelse(votes[i] == m, i - 0x01, 0x00)
        end
        result[site] = reinterpret(eltype(A), merged)
    end
    return result
end
