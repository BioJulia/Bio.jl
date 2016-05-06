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
#
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
#  index :           j - 1              j               j + 1
#   data : ....|xxxxx...........|xxXxxxxxxxxxxxxx|............xxxx|....
# offset :                          |<--- r ----|
#  width :      |<---- 64 ---->| |<---- 64 ---->| |<---- 64 ---->|
#
#  * '.' : unused (4 bits/char)
#  * 'x' : used
#  * 'X' : used and pointed by (j, r)


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

macro dna_str(seq, flags...)
    return DNASequence(remove_newlines(seq))
end

macro rna_str(seq, flags...)
    return RNASequence(remove_newlines(seq))
end

macro aa_str(seq, flags...)
    return AminoAcidSequence(remove_newlines(seq))
end

macro char_str(seq, flags...)
    return CharSequence(remove_newlines(seq))
end


# Constructors
# ------------

seq_data_len{A}(::Type{A}, len::Integer) = cld(len, div(64, bitsof(A)))

function Base.call{A<:Alphabet}(::Type{BioSequence{A}}, len::Integer)
    return BioSequence{A}(Vector{UInt64}(seq_data_len(A, len)), 1:len, false)
end

BioSequence(::Type{DNANucleotide}) = DNASequence()
BioSequence(::Type{RNANucleotide}) = RNASequence()
BioSequence(::Type{AminoAcid}) = AminoAcidSequence()
BioSequence(::Type{Char}) = CharSequence()

function Base.call{A<:Alphabet}(::Type{BioSequence{A}},
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

function Base.call{A}(::Type{BioSequence{A}},
                      other::BioSequence{A},
                      part::UnitRange)
    return BioSequence(other, part)
end

# concatenate chunks
function Base.call{A}(::Type{BioSequence{A}}, chunks::BioSequence{A}...)
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

function Base.(:*){A}(chunk::BioSequence{A}, chunks::BioSequence{A}...)
    return BioSequence{A}(chunk, chunks...)
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

Base.(:^)(chunk::BioSequence, n::Integer) = repeat(chunk, n)

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


# Printers
# --------

Base.summary{A<:DNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "DNA Sequence")
Base.summary{A<:RNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "RNA Sequence")
Base.summary(seq::AminoAcidSequence) = string(length(seq), "aa ", "Amino Acid Sequence")
Base.summary(seq::CharSequence) = string(length(seq), "char ", "Char Sequence")

# pretting printing of sequences
function Base.show(io::IO, seq::BioSequence)
    println(io, summary(seq), ':')
    showcompact(io, seq)
end

function Base.showcompact(io::IO, seq::BioSequence)
    # don't show more than this many characters
    # to avoid filling the screen with junk
    width = displaysize()[2]
    if length(seq) > width
        half = div(width, 2)
        for i in 1:half-1
            print(io, seq[i])
        end
        print(io, '…')
        for i in endof(seq)-half+2:endof(seq)
            print(io, seq[i])
        end
    else
        for x in seq
            write(io, Char(x))
        end
    end
end

# simple printing of sequences
function Base.print(io::IO, seq::BioSequence, width::Integer=50)
    col = 0
    for x in seq
        write(io, Char(x))
        col += 1
        if col == width
            write(io, '\n')
            col = 0
        end
    end
end


# BitIndex
# --------

# useful for accessing bits of sequence's data
immutable BitIndex
    val::Int64
end

index(i::BitIndex) = (i.val >> 6) + 1
offset(i::BitIndex) = i.val & 0b111111
Base.(:+)(i::BitIndex, n::Int) = BitIndex(i.val + n)
Base.(:-)(i::BitIndex, n::Int) = BitIndex(i.val - n)
Base.(:-)(i1::BitIndex, i2::BitIndex) = i1.val - i2.val
Base.(:(==))(i1::BitIndex, i2::BitIndex) = i1.val == i2.val
Base.isless(i1::BitIndex, i2::BitIndex) = isless(i1.val, i2.val)
Base.show(io::IO, i::BitIndex) = print(io, '(', index(i), ", ", offset(i), ')')
Base.start(i::BitIndex) = 1
Base.done(i::BitIndex, s) = s > 2
Base.next(i::BitIndex, s) = ifelse(s == 1, (index(i), 2), (offset(i), 3))

# assumes `i` is positive and `bitsof(A)` is a power of 2
@inline function bitindex{A}(seq::BioSequence{A}, i::Integer)
    return BitIndex((Int(i) + seq.part.start - 2) << trailing_zeros(bitsof(A)))
end


# Basic Operators
# ---------------

Base.length(seq::BioSequence) = length(seq.part)
Base.endof(seq::BioSequence) = length(seq)
Base.isempty(seq::BioSequence) = length(seq) == 0
Base.eltype{A}(::Type{BioSequence{A}}) = eltype(A)

@inline function Base.checkbounds(seq::BioSequence, i::Integer)
    if 1 ≤ i ≤ endof(seq)
        return true
    end
    throw(BoundsError(seq, i))
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

@inline function unsafe_getindex{A}(seq::BioSequence{A}, i::Integer)
    j = bitindex(seq, i)
    @inbounds return decode(A, (seq.data[index(j)] >> offset(j)) & mask(A))
end

function Base.getindex{A}(seq::BioSequence{A}, i::Integer)
    checkbounds(seq, i)
    return unsafe_getindex(seq, i)
end

Base.getindex(seq::BioSequence, part::UnitRange) = BioSequence(seq, part)
Base.sub(seq::BioSequence, part::UnitRange) = BioSequence(seq, part)

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
                         # these types support random access
                         src::Union{AbstractVector,
                                    ASCIIString,
                                    UTF32String,
                                    SubString{ASCIIString},
                                    SubString{UTF32String}},
                         soff::Integer,
                         len::Integer)
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
            x |= enc64(dst, convert(Char, src[i])) << offset(next)
            i += 1
            next += bitsof(A)
        end
        dst.data[j] = x
    end
    return dst
end

function encode_copy!{A}(dst::BioSequence{A},
                         doff::Integer,
                         src::AbstractString,
                         soff::Integer,
                         len::Integer)
    checkbounds(dst, doff:doff+len-1)
    if length(src) < soff + len - 1
        throw(ArgumentError("source string does not contain $len elements from $soff"))
    end

    orphan!(dst)
    next = bitindex(dst, doff)
    stop = bitindex(dst, doff + len)

    # cast first `soff - 1` elements by iteration because strings
    # whose encoding is variable-length don't support random access
    s = start(src)
    for i in 1:soff-1
        _, s = Base.next(src, s)
    end

    while next < stop
        x = UInt64(0)
        j = index(next)
        while index(next) == j && next < stop
            @assert !done(src, s)
            c, s = Base.next(src, s)
            x |= enc64(dst, c) << offset(next)
            next += bitsof(A)
        end
        dst.data[j] = x
    end
    return dst
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

enc64{A}(::BioSequence{A}, x) = UInt64(encode(A, convert(eltype(A), x)))

# Replace a BioSequence's data with a copy, copying only what's needed.
# The user should never need to call this, as it has no outward effect on the
# sequence.
function orphan!{A}(seq::BioSequence{A}, size::Integer=length(seq), force::Bool=false)
    if !seq.shared && !force
        return seq
    end

    j, r = bitindex(seq, 1)
    data = Vector{UInt64}(seq_data_len(A, size))

    if !isempty(seq)
        x = seq.data[j] >> r
        l = min(endof(data), index(bitindex(seq, endof(seq))) - j + 1)
        @inbounds @simd for i in 1:l-1
            y = seq.data[j+i]
            data[i] = x | y << (64 - r)
            x = y >> r
        end
        data[l] = x
    end

    seq.data = data
    seq.part = 1:length(seq)
    seq.shared = false
    return seq
end

function Base.similar{A}(seq::BioSequence{A}, len::Integer=length(seq))
    return BioSequence{A}(len)
end

# actually, users don't need to create a copy of a sequence.
function Base.copy{A}(seq::BioSequence{A})
    # NOTE: no need to set `seq.shared = true` here
    # since `newseq` will be `orphan!`ed soon.
    newseq = BioSequence{A}(seq.data, 1:endof(seq), true)
    orphan!(newseq, length(seq), true)  # force orphan!
    @assert newseq.data !== seq.data
    return newseq
end

# iterator
Base.start(seq::BioSequence) = 1
Base.done(seq::BioSequence, i) = i > endof(seq)
Base.next(seq::BioSequence, i) = unsafe_getindex(seq, i), i + 1

Base.eachindex(seq::BioSequence) = 1:endof(seq)

function Base.(:(==)){A1,A2}(s1::BioSequence{A1}, s2::BioSequence{A2})
    if s1.data === s2.data && s1.part === s2.part && A1 === A2
        return true
    elseif length(s1) != length(s2)
        return false
    end
    for (x, y) in zip(s1, s2)
        if x != y
            return false
        end
    end
    return true
end

function Base.cmp{A1,A2}(s1::BioSequence{A1}, s2::BioSequence{A2})
    for i in 1:min(endof(s1), endof(s2))
        c = cmp(s1[i], s2[i])
        if c != 0
            return c
        end
    end
    return cmp(length(s1), length(s2))
end

Base.isless{A1,A2}(s1::BioSequence{A1}, s2::BioSequence{A2}) = cmp(s1, s2) < 0


# Transformations
# ---------------

function Base.reverse!(seq::BioSequence)
    orphan!(seq)
    for i in 1:div(endof(seq), 2)
        x = unsafe_getindex(seq, i)
        unsafe_setindex!(seq, unsafe_getindex(seq, endof(seq) - i + 1), i)
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
        next = bitindex(seq, endof(seq))
        stop = bitindex(seq, 0)
        i = 0
        while next - stop > 0
            r = offset(next) + $n
            x = seq.data[index(next)] << (64 - r)
            next -= r
            if next - stop > 0
                x |= seq.data[index(next)] >> r
                next -= 64 - r
            end
            data[i+=1] = $nucrev(x)
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

Ambiguous nucleotides are left as-is.
"""
@generated function complement!{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    n = bitsof(A)
    if n == 2
        nuccomp = :nuccomp2
    elseif n == 4
        nuccomp = :nuccomp4
    else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        orphan!(seq)
        next = bitindex(seq, 1)
        stop = bitindex(seq, endof(seq) + 1)
        @inbounds while next < stop
            seq.data[index(next)] = $nuccomp(seq.data[index(next)])
            next += 64
        end
        return seq
    end
end

nuccomp2(x::UInt64) = ~x

@inline function nuccomp4(x::UInt64)
    if x & 0xCCCCCCCCCCCCCCCC != 0
        # ignore ambiguous nucleotides
        m = (x & 0x8888888888888888) | (x & 0x4444444444444444) << 1
        m |= m >> 1 | m >> 2 | m >> 3
        m = ~m
        return (~x & 0x3333333333333333 & m) | (x & ~m)
    else
        return  ~x & 0x3333333333333333
    end
end

"""
    complement(seq)

Make a complement sequence of `seq`.

Ambiguous nucleotides are left as-is.
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
    return reverse_complement!(copy(seq))
end


# Finders
# -------

function Base.findnext(seq::BioSequence, val, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in Int(start):endof(seq)
        x = unsafe_getindex(seq, i)
        if x == v
            return i
        end
    end
    return 0
end

function Base.findprev(seq::BioSequence, val, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in Int(start):-1:1
        x = unsafe_getindex(seq, i)
        if x == v
            return i
        end
    end
    return 0
end

Base.findfirst(seq::BioSequence, val) = findnext(seq, val, 1)
Base.findlast(seq::BioSequence, val)  = findprev(seq, val, endof(seq))

immutable AmbiguousNucleotideIterator{A<:Union{DNAAlphabet,RNAAlphabet}}
    seq::BioSequence{A}
end

ambiguous_positions(seq::BioSequence) = AmbiguousNucleotideIterator(seq)

Base.start(it::AmbiguousNucleotideIterator) = find_next_ambiguous(it.seq, 1)
Base.done(it::AmbiguousNucleotideIterator, nextpos) = nextpos == 0
function Base.next(it::AmbiguousNucleotideIterator, nextpos)
    return nextpos, find_next_ambiguous(it.seq, nextpos + 1)
end

function find_next_ambiguous{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(
        seq::BioSequence{A}, i::Integer)
    # no ambiguity
    return 0
end

function find_next_ambiguous{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(
        seq::BioSequence{A}, i::Integer)
    if i > endof(seq)
        return 0
    end

    # extract ambiguity bits using following bit masks:
    #   * 0x44... = 0b01000100...
    #   * 0x88... = 0b10001000...
    #   * 0xCC... = 0b11001100...
    next = bitindex(seq, max(i, 1))
    stop = bitindex(seq, endof(seq) + 1)
    while next < stop
        m = mask(min(stop - next, 64)) << offset(next)
        x = seq.data[index(next)] & m
        if x & 0xCCCCCCCCCCCCCCCC != 0
            x = (x & 0x8888888888888888) | (x & 0x4444444444444444) << 1
            return (index(next) - 1) << 4 + (trailing_zeros(x) + 1) >> 2
        end
        next += 64 - offset(next)
    end

    # found no ambiguity
    return 0
end


# Mismatch counting
# -----------------

"""
`mismatches(a::NucleotideSequence, b::NucleotideSequence, [nmatches=false])`

Return the number of mismatches between `a` and `b`.

If `a` and `b` are of differing lengths, only the first `min(length(a), length(b))`
nucleotides are compared.

### Arguments
  * `a`: first sequence to compare
  * `b`: second sequence to compare
  * `nmatches`: if true, N matches anything, if false, N matches only itself (false)

### Returns
The number of mismatches
"""
function mismatches end

function mismatches{A<:DNAAlphabet}(a::BioSequence{A},
                                    b::BioSequence{A},
                                    nmatches::Bool=false)
    if nmatches
        return count_mismatches(a, b, nmatches, DNA_N)
    else
        return count_mismatches(a, b)
    end
end

function mismatches{A<:RNAAlphabet}(a::BioSequence{A},
                                    b::BioSequence{A},
                                    nmatches::Bool=false)
    if nmatches
        return count_mismatches(a, b, nmatches, RNA_N)
    else
        return count_mismatches(a, b)
    end
end

function mismatches(a::AminoAcidSequence, b::AminoAcidSequence, xmatches::Bool=false)
    return count_mismatches(a, b, xmatches, AA_X)
end

function count_mismatches(a::BioSequence, b::BioSequence, anychar_matches::Bool, anychar)
    count = 0
    for (x, y) in zip(a, b)
        if anychar_matches
            mismatch = x != y && x != anychar && y != anychar
        else
            mismatch = x != y
        end
        count += mismatch
    end
    return count
end

# fast and exact counting algorithm
@generated function count_mismatches{A<:Union{DNAAlphabet,RNAAlphabet}}(a::BioSequence{A},
                                                                        b::BioSequence{A})
    n = bitsof(A)
    if n == 2
        nucmismatches = :nuc2mismatches
    elseif n == 4
        nucmismatches = :nuc4mismatches
    else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        if length(a) > length(b)
            return count_mismatches(b, a)
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
            mismatches += $nucmismatches(x & m, y & m)
            nexta += k
            nextb += k
        end
        @assert offset(nexta) == 0

        if offset(nextb) == 0  # data are aligned with each other
            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                mismatches += $nucmismatches(x, y)
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                x = a.data[index(nexta)]
                y = b.data[index(nextb)]
                m = mask(stopa - nexta)
                mismatches += $nucmismatches(x & m, y & m)
            end
        elseif nexta < stopa
            y = b.data[index(nextb)]
            nextb += 64

            while stopa - nexta ≥ 64
                x = a.data[index(nexta)]
                z = b.data[index(nextb)]
                y = y >> offset(nextb) | z << (64 - offset(nextb))
                mismatches += $nucmismatches(x, y)
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
                mismatches += $nucmismatches(x & m, y & m)
            end
        end

        return mismatches
    end
end

# mismatch count between two kmers (4 bits)
function nuc4mismatches(x::UInt64, y::UInt64)
    xyxor = x $ y
    mismatches = UInt64(0)
    mismatches |=  xyxor & 0x1111111111111111
    mismatches |= (xyxor & 0x2222222222222222) >> 1
    mismatches |= (xyxor & 0x4444444444444444) >> 2
    mismatches |= (xyxor & 0x8888888888888888) >> 3
    return count_ones(mismatches)
end

# mismatch count between two kmers (2 bits)
function nuc2mismatches(x::UInt64, y::UInt64)
    xyxor = x $ y
    mismatches = UInt64(0)
    mismatches |=  xyxor & 0x5555555555555555
    mismatches |= (xyxor & 0xAAAAAAAAAAAAAAAA) >> 1
    return count_ones(mismatches)
end
