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
#     j, r = bitsid(seq, i)
#     decode(A, UInt8((seq.data[j] >> r) & mask(A)))
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
end

# type aliases
typealias DNASequence       BioSequence{DNAAlphabet{4}}
typealias RNASequence       BioSequence{RNAAlphabet{4}}
typealias AminoAcidSequence BioSequence{AminoAcidAlphabet}


# Constructors
# ------------

# This is the body of the BioSequence constructor below. It's separated into a
# macro so we can generate two versions depending on wether the `unsafe` flag is
# set. In this macro, `srcdata[startpos:stoppos]` would be encoded into
# `seqdata` as alphabet `A`.
macro encode_seq(nt_convert_expr)
    quote
        j = startpos
        @inbounds begin
            hasinvalid = false
            for i in 1:endof(seqdata)
                n = bitsof(A)
                shift = 0
                data_i = UInt64(0)
                while shift < 64 && j <= stoppos
                    c = srcdata[j]
                    bioc = $(nt_convert_expr)
                    hasinvalid |= !isvalid(bioc)
                    data_i |= UInt64(encode(A, bioc)) << shift
                    j += 1
                    shift += n
                end
                seqdata[i] = data_i
            end

            if hasinvalid
                # figure out what the first bad character was.
                for i in startpos:stoppos
                    c = srcdata[j]
                    bioc = $(nt_convert_expr)
                    if !isvalid(bioc)
                        error(c, " is not a valid character in ", A)
                    end
                end
            end
        end
    end
end

seq_data_len{A}(::Type{A}, len::Integer) = cld(len, div(64, bitsof(A)))

function Base.call{A<:Alphabet}(::Type{BioSequence{A}})
    return BioSequence{A}(UInt64[], 1:0, false)
end

function Base.call{A<:Alphabet}(::Type{BioSequence{A}}, len::Integer)
    return BioSequence{A}(Vector{UInt64}(seq_data_len(A, len)), 1:len, false)
end

BioSequence(::Type{DNANucleotide}) = DNASequence()
BioSequence(::Type{RNANucleotide}) = RNASequence()
BioSequence(::Type{AminoAcid}) = AminoAcidSequence()

function Base.call{A<:Alphabet}(::Type{BioSequence{A}},
                                seq::Union{AbstractString,Vector{UInt8}},
                                startpos::Integer=1,
                                stoppos::Integer=endof(seq);
                                unsafe::Bool=false)
    len = stoppos - startpos + 1
    seqdata = zeros(UInt64, seq_data_len(A, len))
    srcdata = seq
    T = eltype(A)
    if unsafe
        @encode_seq unsafe_ascii_byte_to_nucleotide(T, c)
    else
        @encode_seq convert(T, Char(c))
    end
    return BioSequence{A}(seqdata, 1:len, false)
end

function Base.call{A<:DNAAlphabet}(::Type{BioSequence{A}},
                                   seq::AbstractVector{DNANucleotide},
                                   startpos::Integer=1,
                                   stoppos::Integer=endof(seq);
                                   unsafe::Bool=false)
    return make_from_vector(A, seq, startpos, stoppos, unsafe)
end

function Base.call{A<:RNAAlphabet}(::Type{BioSequence{A}},
                                   seq::AbstractVector{RNANucleotide},
                                   startpos::Integer=1,
                                   stoppos::Integer=endof(seq);
                                   unsafe::Bool=false)
    return make_from_vector(A, seq, startpos, stoppos, unsafe)
end

function Base.call(::Type{AminoAcidSequence},
                   seq::AbstractVector{AminoAcid},
                   startpos::Integer=1,
                   stoppos::Integer=endof(seq);
                   unsafe::Bool=false)
    return make_from_vector(AminoAcidAlphabet, seq, startpos, stoppos, unsafe)
end

function make_from_vector{A,T}(::Type{A}, srcdata::AbstractVector{T}, startpos, stoppos, unsafe)
    len = stoppos - startpos + 1
    seqdata = zeros(UInt64, seq_data_len(A, len))
    if unsafe
        @encode_seq c
    else
        @encode_seq begin
            if !isvalid(c)
                throw(ArgumentError("invalid $T"))
            end
            c
        end
    end
    return BioSequence{A}(seqdata, 1:len, false)
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

function BioSequence{A}(chunks::BioSequence{A}...)
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

# conversion between DNA and RNA
for (A1, A2) in [(DNAAlphabet, RNAAlphabet), (RNAAlphabet, DNAAlphabet)], n in (2, 4)
    # NOTE: assumes that binary representation is identical between DNA and RNA
    @eval function Base.convert(::Type{BioSequence{$(A1{n})}}, seq::BioSequence{$(A2{n})})
        newseq = BioSequence{$(A1{n})}(seq.data, seq.part, true)
        seq.shared = true
        return newseq
    end
end

# from a vector
Base.convert{A<:DNAAlphabet}(::Type{BioSequence{A}}, seq::AbstractVector{DNANucleotide}) =
    BioSequence{A}(seq, 1, endof(seq))
Base.convert{A<:RNAAlphabet}(::Type{BioSequence{A}}, seq::AbstractVector{RNANucleotide}) =
    BioSequence{A}(seq, 1, endof(seq))
Base.convert(::Type{AminoAcidSequence}, seq::AbstractVector{AminoAcid}) =
    AminoAcidSequence(seq, 1, endof(seq))

# to a vector
Base.convert(::Type{Vector}, seq::BioSequence) = collect(seq)
Base.convert{A<:DNAAlphabet}(::Type{Vector{DNANucleotide}}, seq::BioSequence{A}) = collect(seq)
Base.convert{A<:RNAAlphabet}(::Type{Vector{RNANucleotide}}, seq::BioSequence{A}) = collect(seq)
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = collect(seq)

# from/to a string
Base.convert{S<:AbstractString}(::Type{S}, seq::BioSequence) = convert(S, [Char(x) for x in seq])
Base.convert{S<:AbstractString,A}(::Type{BioSequence{A}}, seq::S) = BioSequence{A}(seq)


# Printers
# --------

Base.summary{A<:DNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "DNA Sequence")
Base.summary{A<:RNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "RNA Sequence")
Base.summary(seq::AminoAcidSequence) = string(length(seq), "aa ", "Amino Acid Sequence")

# pretting printing of sequences
function Base.show(io::IO, seq::BioSequence)
    println(io, summary(seq), ':')
    showcompact(io, seq)
end

function Base.showcompact(io::IO, seq::BioSequence)
    # don't show more than this many characters
    # to avoid filling the screen with junk
    width = Base.tty_size()[2] - 2
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


# Basic Operators
# ---------------

Base.length(seq::BioSequence) = length(seq.part)
Base.endof(seq::BioSequence) = length(seq)
Base.isempty(seq::BioSequence) = length(seq) == 0
Base.eltype{A}(seq::BioSequence{A}) = eltype(A)

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

# Returns `(<element index>, <bits offset>)`.
# NOTE: This function assumes `i ≥ 1`.
@inline function bitsid{A}(seq::BioSequence{A}, i::Integer)
    n = bitsof(A)
    d, r = divrem(i + seq.part.start - 2, div(64, n))
    return d + 1, r * n
end

@inline function bitsid{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(seq::BioSequence{A},
                                                                 i::Integer)
    d, r = divrem32(i + seq.part.start - 2)
    return d + 1, 2r
end

@inline function bitsid{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(seq::BioSequence{A},
                                                                 i::Integer)
    d, r = divrem16(i + seq.part.start - 2)
    return d + 1, 4r
end

@inline function bitsid(seq::AminoAcidSequence, i::Integer)
    d, r = divrem8(i + seq.part.start - 2)
    return d + 1, 8r
end

#=
# Unfortunately, combination of `@inline` and `@generated` won't work.
@inline @generated function bitsid{A}(seq::BioSequence{A}, i::Integer)
    n = bitsof(A)
    index = :(Int(i) + seq.part.start - 2)

    if n == 2
        divrem = :(divrem32($index))
    elseif n == 4
        divrem = :(divrem16($index))
    elseif n == 8
        divrem = :( divrem8($index))
    else
        divrem = :(  divrem($index, $(div(64, n))))
    end

    quote
        d, r = $divrem
        # (element index, bits offset)
        return d + 1, r * $n
    end
end
=#

mask(n::Integer) = (UInt64(1) << n) - 1
mask{A<:Alphabet}(::Type{A}) = mask(bitsof(A))

divrem8(i::Int)  = i >> 3, i & 0b111
divrem16(i::Int) = i >> 4, i & 0b1111
divrem32(i::Int) = i >> 5, i & 0b11111

@inline function unsafe_getindex{A}(seq::BioSequence{A}, i::Integer)
    j, r = bitsid(seq, i)
    @inbounds return decode(A, UInt8((seq.data[j] >> r) & mask(A)))
end

function Base.getindex{A}(seq::BioSequence{A}, i::Integer)
    checkbounds(seq, i)
    return unsafe_getindex(seq, i)
end

Base.getindex(seq::BioSequence, part::UnitRange) = BioSequence(seq, part)
Base.sub(seq::BioSequence, part::UnitRange)      = BioSequence(seq, part)

function Base.repeat{A}(chunk::BioSequence{A}, n::Integer)
    seq = BioSequence{A}(length(chunk) * n)
    offset = 1
    for _ in 1:n
        copy!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

Base.(:*){A}(chunk1::BioSequence{A}, chunks::BioSequence{A}...) =
    BioSequence(chunk1, chunks...)

Base.(:^)(chunk::BioSequence, n::Integer) = repeat(chunk, n)

function Base.setindex!(seq::BioSequence, x, i::Integer)
    checkbounds(seq, i)
    orphan!(seq)
    return unsafe_setindex!(seq, x, i)
end

@generated function Base.setindex!{A,T<:Integer}(seq::BioSequence{A},
                                                 x,
                                                 locs::AbstractVector{T})
    if locs <: AbstractVector{Bool}
        quote
            checkbounds(seq, locs)
            bin = enc(seq, x)
            orphan!(seq)
            i = j = 0
            while (i = findnext(locs, i + 1)) > 0
                encoded_setindex!(seq, bin, i)
            end
            return seq
        end
    else
        quote
            checkbounds(seq, locs)
            bin = enc(seq, x)
            orphan!(seq)
            for i in locs
                encoded_setindex!(seq, bin, i)
            end
            return seq
        end
    end
end

Base.setindex!{A}(seq::BioSequence{A}, x, ::Colon) = setindex!(seq, x, 1:endof(seq))

@generated function Base.setindex!{A,T<:Integer}(seq::BioSequence{A},
                                                 other::BioSequence{A},
                                                 locs::AbstractVector{T})
    if locs <: UnitRange
        quote
            checkbounds(seq, locs)
            checkdimension(other, locs)
            return copy!(seq, locs.start, other, 1)
        end
    elseif locs <: AbstractVector{Bool}
        quote
            checkbounds(seq, locs)
            checkdimension(other, locs)
            orphan!(seq)
            i = j = 0
            while (i = findnext(locs, i + 1)) > 0
                unsafe_setindex!(seq, other[j+=1], i)
            end
            return seq
        end
    else
        quote
            checkbounds(seq, locs)
            checkdimension(other, locs)
            orphan!(seq)
            for (i, x) in zip(locs, other)
                unsafe_setindex!(seq, x, i)
            end
            return seq
        end
    end
end

function Base.setindex!{A}(seq::BioSequence{A}, other::BioSequence{A}, ::Colon)
    return setindex!(seq, other, 1:endof(seq))
end

@inline function unsafe_setindex!{A}(seq::BioSequence{A}, x, i::Integer)
    bin = enc(seq, x)
    return encoded_setindex!(seq, UInt64(bin), i)
end

@inline function encoded_setindex!{A}(seq::BioSequence{A}, bin::UInt64, i::Integer)
    j, r = bitsid(seq, i)
    m = mask(A)
    @inbounds seq.data[j] = (bin << r) | (seq.data[j] & ~(m << r))
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
    bin = enc(seq, x)
    resize!(seq, length(seq) + 1)
    encoded_setindex!(seq, bin, endof(seq))
    return seq
end

function Base.unshift!{A}(seq::BioSequence{A}, x)
    bin = enc(seq, x)
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
    bin = enc(seq, x)
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

function Base.copy!{A}(seq::BioSequence{A},
                       srcdata::Vector{UInt8},
                       startpos::Integer, stoppos::Integer)
    n = stoppos - startpos + 1
    len = seq_data_len(A, n)
    seqdata = seq.data
    if length(seqdata) < len
        resize!(seq.data, len)
    end
    fill!(seq.data, 0)
    seq.part = 1:n
    @encode_seq convert(eltype(A), Char(c))
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

    j1, r1 = bitsid(dst, doff)
    j2, r2 = bitsid(src, soff)
    rest = len * bitsof(A)

    while rest > 0
        # move `k` bits from `src` to `dst`
        x = dst.data[j1]
        y = src.data[j2]
        if r1 < r2
            y >>= r2 - r1
            k = min(64 - r2, rest)
        else
            y <<= r1 - r2
            k = min(64 - r1, rest)
        end
        m = mask(k) << r1
        dst.data[j1] = y & m | x & ~m

        if r1 + k ≥ 64
            j1 += 1
            r1 = 0
        else
            r1 += k
        end
        if r2 + k ≥ 64
            j2 += 1
            r2 = 0
        else
            r2 += k
        end
        rest -= k
    end

    return dst
end

@inline function enc{A}(seq::BioSequence{A}, x)
    return UInt64(encode(A, convert(eltype(A), x)))
end

# Replace a BioSequence's data with a copy, copying only what's needed.
# The user should never need to call this, as it has no outward effect on the
# sequence.
function orphan!{A}(seq::BioSequence{A}, size::Integer=length(seq), force::Bool=false)
    if !seq.shared && !force
        return seq
    end

    j, r = bitsid(seq, 1)
    data = Vector{UInt64}(seq_data_len(A, size))

    if !isempty(seq)
        x = seq.data[j] >> r
        l = min(endof(data), bitsid(seq, endof(seq))[1] - j + 1)
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

# actually, users don't need to create a copy of a sequence.
function Base.copy(seq::BioSequence)
    newseq = seq[1:end]
    orphan!(newseq)
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

function Base.reverse!{A}(seq::BioSequence{A})
    orphan!(seq)
    for i in 1:div(endof(seq), 2)
        seq[i], seq[end-i+1] = seq[end-i+1], seq[i]
    end
    return seq
end

function Base.reverse{A}(seq::BioSequence{A})
    return reverse!(copy(seq))
end

@generated function Base.reverse{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    n = bitsof(A)
    if n == 2
        nucrev = :nucrev2
    elseif n == 4
        nucrev = :nucrev4
    else
        error("n ∉ (2, 4)")
    end

    quote
        data = zeros(UInt64, seq_data_len(A, length(seq)))
        j, r = bitsid(seq, endof(seq) + 1)
        if r == 0
            @inbounds for j′ in 1:endof(data)
                j -= 1
                data[j′] = $nucrev(seq.data[j])
            end
        else
            x = (seq.data[j] & mask(r)) << (64 - r)
            @inbounds for j′ in 1:endof(data)-1
                j -= 1
                x′ = seq.data[j]
                data[j′] = $nucrev(x | (x′ >> r))
                x = x′ << (64 - r)
            end
            data[end] = $nucrev(x)
        end
        return BioSequence{A}(data, 1:length(seq), false)
    end
end

function Base.complement!{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(seq::BioSequence{A})
    if isempty(seq)
        return seq
    end
    orphan!(seq)
    j1, r1 = bitsid(seq, 1)
    j2, r2 = bitsid(seq, endof(seq) + 1)
    data_j = seq.data[j1]
    if j1 == j2
        m = mask(r2 - r1) << r1
        seq.data[j1] = (~data_j & m) | (data_j & ~m)
    else
        m = mask(64 - r1) << r1
        seq.data[j1] = (~data_j & m) | (data_j & ~m)
        for j in j1+1:j2-1
            seq.data[j] = ~seq.data[j]
        end
        if r2 > 0
            data_j = seq.data[j2]
            m = mask(r2 + 2)
            seq.data[j2] = (~data_j & m) | (data_j & ~m)
        end
    end
    return seq
end

# Compute complement of `x`'s bits selected by `mask`.
# NOTE: 4-bit encoding is assumed and ambiguous nucleotides are left untouched.
macro complement(x, mask)
    quote
        if hasambiguous($x)
            m = ($x & 0x8888888888888888) | (($x & 0x4444444444444444) << 1)
            m |= m >> 1 | m >> 2 | m >> 3
            m = ~m & $mask
            (~$x & 0x3333333333333333 & m) | ($x & ~m)
        else
            m = $mask
            (~$x & 0x3333333333333333 & m) | ($x & ~m)
        end
    end
end

function Base.complement!(seq::Union{DNASequence,RNASequence})
    if isempty(seq)
        return seq
    end
    orphan!(seq)
    j1, r1 = bitsid(seq, 1)
    j2, r2 = bitsid(seq, endof(seq) + 1)

    # special case: whole sequence is within an element (i.e. `seq.data[j1]`)
    if j1 == j2
        seq.data[j1] = @complement seq.data[j1] (mask(r2) - mask(r1))
        return seq
    end

    seq.data[j1] = @complement seq.data[j1] ~mask(r1)
    @inbounds for j in j1+1:j2-1
        seq.data[j] = @complement seq.data[j] 0xFFFFFFFFFFFFFFFF
    end
    if r2 > 0
        seq.data[j2] = @complement seq.data[j2] mask(r2)
    end
    return seq
end

hasambiguous(x::UInt64) = x & 0xCCCCCCCCCCCCCCCC != 0

function Base.complement{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(copy(seq))
end

function reverse_complement!{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(reverse!(seq))
end

function reverse_complement{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return reverse_complement!(copy(seq))
end

function nucrev2(x::UInt64)
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>>  2
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >>> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >>> 32
     return x
end

function nucrev4(x::UInt64)
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >>> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >>> 32
     return x
end

immutable AmbiguousNucleotideIterator{A<:Union{DNAAlphabet,RNAAlphabet}}
    seq::BioSequence{A}
end

function npositions{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return AmbiguousNucleotideIterator(seq)
end

Base.start(it::AmbiguousNucleotideIterator) = find_next_ambiguous(it.seq, 1)
Base.done(it::AmbiguousNucleotideIterator, nextpos) = nextpos == 0
function Base.next(it::AmbiguousNucleotideIterator, nextpos)
    return nextpos, find_next_ambiguous(it.seq, nextpos + 1)
end

function find_next_ambiguous{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(seq::BioSequence{A}, i::Integer)
    # no ambiguity
    return 0
end

# Note: `i` is absolute index
function find_next_ambiguous{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(seq::BioSequence{A}, i::Integer)
    if i > endof(seq)
        return 0
    end
    j, r = divrem16(Int(max(i, 1) + seq.part.start - 2))
    j += 1
    data_j = seq.data[j] & ~mask(4r)
    # extract ambiguity bits using following bit masks:
    #   * 0x44... = 0b01000100...
    #   * 0x88... = 0b10001000...
    #   * 0xCC... = 0b11001100...
    amb = data_j & 0xCCCCCCCCCCCCCCCC
    if amb != 0
        amb = ((amb & 0x8888888888888888) >> 1) | (amb & 0x4444444444444444)
        return ((j - 1) << 4) + ((trailing_zeros(amb) - 2) >> 2) + 1
    end
    for j in j+1:endof(seq.data)
        data_j = seq.data[j]
        amb = data_j & 0xCCCCCCCCCCCCCCCC
        if amb != 0
            amb = ((amb & 0x8888888888888888) >> 1) | (amb & 0x4444444444444444)
            return ((j - 1) << 4) + ((trailing_zeros(amb) - 2) >> 2) + 1
        end
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
        error("n ∉ (2, 4)")
    end

    quote
        if length(a) > length(b)
            return count_mismatches(b, a)
        end

        ja, ra = bitsid(a, 1)
        j2, r2 = bitsid(a, endof(a) + 1)
        jb, rb = bitsid(b, 1)
        mismatches = 0

        m = mask($n)
        while ra != 0 && (ja, ra) < (j2, r2)
            x = (a.data[ja] >> ra) & m
            y = (b.data[jb] >> rb) & m
            mismatches += x != y
            if (ra += $n) ≥ 64
                ja += 1
                ra = 0
            end
            if (rb += $n) ≥ 64
                jb += 1
                rb = 0
            end
        end

        if ja == j2 && ra == r2
            return mismatches
        end

        if rb == 0
            # fortunately, data are aligned with each other
            while ja < j2
                x = a.data[ja]
                y = b.data[jb]
                mismatches += $nucmismatches(x, y)
                ja += 1
                jb += 1
            end
            if r2 > 0
                x = a.data[ja]
                y = b.data[jb]
                m = mask(r2)
                mismatches += $nucmismatches(x & m, y & m)
            end
        else
            y = b.data[jb] >> rb
            m = mask(64 - rb)
            while ja < j2
                jb += 1
                x = a.data[ja]
                y′ = b.data[jb]
                mismatches += $nucmismatches(x, y & m | (y′ << (64 - rb)) & ~m)
                y = y′ >> rb
                ja += 1
            end

            m = mask($n)
            while (ja, ra) < (j2, r2)
                x = (a.data[ja] >> ra) & m
                y = (b.data[jb] >> rb) & m
                mismatches += x != y
                if (ra += $n) ≥ 64
                    ja += 1
                    ra = 0
                end
                if (rb += $n) ≥ 64
                    jb += 1
                    rb = 0
                end
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
