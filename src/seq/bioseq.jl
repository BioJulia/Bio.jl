# BioSequence
# ===========
#
# A general purpose biological sequence representation.
#
# TODO: Commens below may not reflect the current changes.
#
# Sequences are either explicitly immutable or mutable. Immutable sequences have
# the benefit that subsequence operations (`seq[i:j]`) are cheap, since they
# share the underlying data and do not copy anything. Mutable sequences may be
# mutated, but as a consequence subsequences copy data rather than reference it.
#
# Sequences can be converted between mutable and immutable using `mutable!` and
# `immutable!`, respectively. Converting from mutable to immutable is cheap: it
# only flips the `mutable` flag. The converse, immutable to mutable, is cheap
# *if* the sequence is not a subsequence and has no subsequences, otherwise it
# must make a copy of the data.
#
#      the behavior of subsequence syntax
#   |           | seq[i:j] | sub(seq, i:j) |
#   |-----------|----------|---------------|
#   | immutable |   view   |     view      |
#   |  mutable  |   copy   |     view      |
#

"""
Biological sequence data structure indexed by an alphabet type `A`.
"""
type BioSequence{A<:Alphabet} <: Sequence
    data::Vector{UInt64}  # encoded sequence
    part::UnitRange{Int}  # interval within `data` defining the (sub)sequence
    mutable::Bool         # true if and only if the sequence can be safely mutated
    shared::Bool          # true if and only if `data` is shared between sequences
end

# type aliases
typealias DNASequence BioSequence{DNAAlphabet{4}}
typealias RNASequence BioSequence{RNAAlphabet{4}}
typealias AminoAcidSequence BioSequence{AminoAcidAlphabet}


# Constructors
# ------------

# This is the body of the BioSequence constructor below. It's separated into a
# macro so we can generate two versions depending on wether the `unsafe` flag is
# set. In this macro, `strdata[startpos:stoppos]` would be encoded into
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
                    c = strdata[j]
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
                    c = strdata[j]
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

function Base.call{A<:Alphabet}(::Type{BioSequence{A}}; mutable::Bool=true)
    return BioSequence{A}(UInt64[], 1:0, mutable, false)
end

function Base.call{A<:Alphabet}(::Type{BioSequence{A}},
                                seq::Union{AbstractString,Vector{UInt8}},
                                startpos::Integer=1,
                                stoppos::Integer=endof(seq);
                                unsafe::Bool=false, mutable::Bool=false)
    len = stoppos - startpos + 1
    seqdata = zeros(UInt64, seq_data_len(A, len))
    strdata = seq
    T = eltype(A)
    if unsafe
        @encode_seq unsafe_ascii_byte_to_nucleotide(T, c)
    else
        @encode_seq convert(T, Char(c))
    end
    return BioSequence{A}(seqdata, 1:len, mutable, false)
end

for A in (DNAAlphabet, RNAAlphabet)
    @eval function Base.call{A<:$(A)}(::Type{BioSequence{A}},
                                    seq::AbstractVector{$(eltype(A))},
                                    startpos::Integer=1,
                                    stoppos::Integer=endof(seq);
                                    unsafe::Bool=false, mutable::Bool=false)
        len = stoppos - startpos + 1
        seqdata = zeros(UInt64, seq_data_len(A, len))
        strdata = seq
        if unsafe
            @encode_seq c
        else
            @encode_seq begin
                if !isvalid(c)
                    error("invalid nucleotide")
                end
                if A == $(A){2} && isambiguous(c)
                    error("cannot store ambiguous nucleotide")
                end
                c
            end
        end
        return BioSequence{A}(seqdata, 1:len, mutable, false)
    end
end

function Base.call(::Type{AminoAcidSequence},
                   seq::AbstractVector{AminoAcid},
                   startpos::Integer=1,
                   stoppos::Integer=endof(seq);
                   unsafe::Bool=false, mutable::Bool=false)
    len = stoppos - startpos + 1
    seqdata = zeros(UInt64, seq_data_len(AminoAcidAlphabet, len))
    strdata = seq
    A = AminoAcidAlphabet
    if unsafe
        @encode_seq c
    else
        @encode_seq begin
            if !isvalid(c)
                error("invalid amino acid")
            end
            c
        end
    end
    return AminoAcidSequence(seqdata, 1:len, mutable, false)
end

function BioSequence{A}(other::BioSequence{A},
                        part::UnitRange;
                        mutable::Bool=other.mutable)
    checkbounds(other, part)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = BioSequence{A}(other.data, start:stop, true, true)
    if mutable || other.mutable
        orphan!(subseq)
    else
        subseq.shared = other.shared = true
    end
    subseq.mutable = mutable
    return subseq
end

function Base.call{A}(::Type{BioSequence{A}},
                      other::BioSequence{A},
                      part::UnitRange;
                      mutable::Bool=other.mutable)
    return BioSequence(other, part, mutable=mutable)
end

function BioSequence{A}(chunks::BioSequence{A}...)
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seqdata = zeros(UInt64, seq_data_len(A, len))
    newseq = BioSequence{A}(seqdata, 1:len, false, false)
    offset = 1
    for chunk in chunks
        unsafe_copy!(newseq, offset, chunk)
        offset += length(chunk)
    end
    return newseq
end

# conversion between DNA and RNA
for (A1, A2) in [(DNAAlphabet, RNAAlphabet), (RNAAlphabet, DNAAlphabet)]
    for n in (2, 4)
        @eval function Base.convert(::Type{BioSequence{$(A1{n})}}, seq::BioSequence{$(A2{n})})
            newseq = BioSequence{$(A1{n})}(seq.data, seq.part, true, true)
            orphan!(newseq)
            newseq.mutable = seq.mutable
            return newseq
        end
    end
end

Base.convert{A<:DNAAlphabet}(::Type{BioSequence{A}}, seq::AbstractVector{DNANucleotide}) =
    BioSequence{A}(seq, 1, endof(seq))
Base.convert{A<:RNAAlphabet}(::Type{BioSequence{A}}, seq::AbstractVector{RNANucleotide}) =
    BioSequence{A}(seq, 1, endof(seq))
Base.convert(::Type{AminoAcidSequence}, seq::AbstractVector{AminoAcid}) =
    AminoAcidSequence(seq, 1, endof(seq))

Base.convert(::Type{Vector}, seq::BioSequence) = collect(seq)
Base.convert{A<:DNAAlphabet}(::Type{Vector{DNANucleotide}}, seq::BioSequence{A}) = collect(seq)
Base.convert{A<:RNAAlphabet}(::Type{Vector{RNANucleotide}}, seq::BioSequence{A}) = collect(seq)
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = collect(seq)

Base.convert{S<:AbstractString}(::Type{S}, seq::BioSequence) = convert(S, [Char(x) for x in seq])
Base.convert{S<:AbstractString,A<:Alphabet}(::Type{BioSequence{A}}, seq::S) = BioSequence{A}(seq)


# Printers
# --------

Base.summary{A<:DNAAlphabet}(seq::BioSequence{A}) =
    string(length(seq), "nt ", seq.mutable ? "Mutable " : "", "DNA Sequence")
Base.summary{A<:RNAAlphabet}(seq::BioSequence{A}) =
    string(length(seq), "nt ", seq.mutable ? "Mutable " : "", "RNA Sequence")
Base.summary(seq::AminoAcidSequence) =
    string(length(seq), "aa ", seq.mutable ? "Mutable " : "", "Amino Acid Sequence")

# pretting printing of sequences
function Base.show(io::IO, seq::BioSequence)
    println(io, summary(seq), ':')
    # don't show more than this many characters to avoid filling the screen
    # with junk
    width = Base.tty_size()[2] - 2
    if length(seq) > width
        half = div(width, 2)
        print(io, seq[1:half-1], -1)
        write(io, " … ")
        print(io, seq[end-half+2:end], -1)
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

@inline function Base.checkbounds(seq::BioSequence, part::UnitRange)
    if 1 ≤ part.start && part.stop ≤ endof(seq)
        return true
    end
    throw(BoundsError(seq, part))
end

@inline function checkmutability(seq::BioSequence)
    if seq.mutable
        return true
    end
    error("attempt to modify immutable sequence")
end

# NOTE: This function assumes `i ≥ 1`.
@inline function bitsid{A}(seq::BioSequence{A}, i::Integer)
    n = bitsof(A)
    d, r = divrem(i + seq.part.start - 2, div(64, n))
    # (element index, bits offset)
    return d + 1, r * n
end

@inline function bitsid{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(seq::BioSequence{A}, i::Integer)
    d, r = divrem32(i + seq.part.start - 2)
    return d + 1, 2r
end

@inline function bitsid{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(seq::BioSequence{A}, i::Integer)
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

# always create a view without copying
function Base.sub{A}(seq::BioSequence{A}, part::UnitRange)
    checkbounds(seq, part)
    start = seq.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = BioSequence{A}(seq.data, start:stop, seq.mutable, true)
    seq.shared = true
    return subseq
end

# util to define a setindex! method
macro define_setindex!(A)
    quote
        @eval function Base.setindex!(seq::BioSequence{$(A)}, x::$(eltype(A)), i::Integer)
            checkmutability(seq)
            checkbounds(seq, i)
            j, r = bitsid(seq, i)
            seq.data[j] = (seq.data[j] & ~(mask($(A)) << r)) | UInt64(encode($(A), x)) << r
            return seq
        end
    end
end

for A in (DNAAlphabet{2}, DNAAlphabet{4},
          RNAAlphabet{2}, RNAAlphabet{4},
          AminoAcidAlphabet)
    @eval begin
        @define_setindex! $(A)
        function Base.setindex!(seq::BioSequence{$(A)}, x::Char, i::Integer)
            return setindex!(seq, $(eltype(A))(x), i)
        end
    end
end

Base.(:*){A}(chunk1::BioSequence{A}, chunks::BioSequence{A}...) =
    BioSequence(chunk1, chunks...)

function Base.repeat{A}(chunk::BioSequence{A}, n::Integer)
    len = length(chunk) * n
    seqdata = zeros(UInt64, seq_data_len(A, len))
    newseq = BioSequence{A}(seqdata, 1:len, false, false)
    offset = 1
    for _ in 1:n
        unsafe_copy!(newseq, offset, chunk)
        offset += length(chunk)
    end
    return newseq
end

Base.(:^)(chunk::BioSequence, n::Integer) = repeat(chunk, n)

function Base.copy!{A}(seq::BioSequence{A},
                       strdata::Vector{UInt8},
                       startpos::Integer, stoppos::Integer)
    checkmutability(seq)
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

function unsafe_copy!{A}(dst::BioSequence{A}, pos::Int, src::BioSequence{A})
    n = bitsof(A)
    j1, r1 = bitsid(dst, pos)
    j2, r2 = bitsid(src, 1)
    rest = length(src) * n
    while rest > 0
        m = mask(min(rest, 64))
        if r1 == r2 == 0
            dst.data[j1] = src.data[j2] & m
        else
            dst.data[j1] |= ((src.data[j2] >> r2) & m) << r1
        end
        if r1 > r2
            k = 64 - r1
            j1 += 1
            r1 = 0
            r2 += k
        elseif r2 > r1
            k = 64 - r2
            j2 += 1
            r2 = 0
            r1 += k
        else # r1 == r2
            k = 64 - r1
            r1 += k
            r2 += k
            if r1 ≥ 64
                r1 = 0
                r2 = 0
                j1 += 1
                j2 += 1
            end
        end
        rest -= k
    end
    return dst
end

# Replace a BioSequence's data with a copy, copying only what's needed.
# The user should never need to call this, as it has no outward effect on the
# sequence, but it makes functions like mismatch easier and faster if can assume
# a sequence is aligned with its data.
function orphan!{A}(seq::BioSequence{A})
    # no need to orphan data from an immutable sequence
    @assert seq.mutable

    data = zeros(UInt64, seq_data_len(A, length(seq)))
    j, shift = divrem(seq.part.start - 1, div(64, bitsof(A)))
    j += 1
    shift *= bitsof(A)

    @inbounds for i in 1:length(data)
        data[i] |= seq.data[j] >> shift
        j += 1
        if j > endof(seq.data)
            break
        end
        data[i] |= seq.data[j] << (64 - shift)
    end

    seq.data = data
    seq.part = 1:length(seq.part)
    seq.shared = false
    return seq
end

function Base.copy(seq::BioSequence)
    if !seq.mutable
        return seq
    end
    return docopy(seq)
end

function docopy{A}(seq::BioSequence{A})
    newseq = BioSequence{A}(seq.data, seq.part, true, false)
    orphan!(newseq)
    @assert newseq.data !== seq.data
    newseq.mutable = seq.mutable
    return newseq
end

# iterator
Base.start(seq::BioSequence) = 1
Base.done(seq::BioSequence, i) = i > endof(seq)
Base.next(seq::BioSequence, i) = unsafe_getindex(seq, i), i + 1

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

ismutable(seq::BioSequence) = seq.mutable

function mutable!(seq::BioSequence)
    if ismutable(seq)
        return seq
    end
    seq.mutable = true
    orphan!(seq)
    return seq
end

function immutable!(seq::BioSequence)
    if !ismutable(seq)
        return seq
    end
    seq.mutable = false
    if seq.shared
        orphan!(seq)
    end
    return seq
end


# Transformations
# ---------------

function Base.reverse!{A}(seq::BioSequence{A})
    checkmutability(seq)
    for i in 1:div(endof(seq), 2)
        seq[i], seq[end-i+1] = seq[end-i+1], seq[i]
    end
    return seq
end

function Base.reverse{A}(seq::BioSequence{A})
    newseq = docopy(seq)
    mutable!(newseq)
    reverse!(newseq)
    newseq.mutable = seq.mutable
    return newseq
end

function Base.complement!{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(seq::BioSequence{A})
    checkmutability(seq)
    if isempty(seq)
        return seq
    end
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
    checkmutability(seq)
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
    newseq = docopy(seq)
    mutable!(newseq)
    complement!(newseq)
    newseq.mutable = seq.mutable
    return newseq
end

function reverse_complement!{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    checkmutability(seq)
    return complement!(reverse!(seq))
end

function reverse_complement{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    newseq = docopy(seq)
    mutable!(newseq)
    complement!(reverse!(newseq))
    newseq.mutable = seq.mutable
    return newseq
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
function mismatches{A<:Union{DNAAlphabet,RNAAlphabet}}(a::BioSequence{A},
                                                       b::BioSequence{A},
                                                       nmatches::Bool=false)
    count = 0
    # TODO: nmatches
    for (x, y) in zip(a, b)
        if x != y
            count += 1
        end
    end
    return count
end

# mismatch count between two kmers
function nucmismatches(x::UInt64, y::UInt64)
    xyxor = x $ y
    return count_ones((xyxor & 0x5555555555555555) | ((xyxor & 0xAAAAAAAAAAAAAAAA) >>> 1))
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
