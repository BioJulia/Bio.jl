
include("bioseq/type.jl")
include("bioseq/constructors.jl")
include("bioseq/conversion.jl")
include("bioseq/stringliterals.jl")
include("bioseq/operators.jl")
include("bioseq/transformations.jl")

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





# Ambiguous nucleotides iterator
# ------------------------------

immutable AmbiguousNucleicAcidIterator{A<:Union{DNAAlphabet,RNAAlphabet}}
    seq::BioSequence{A}
end

ambiguous_positions(seq::BioSequence) = AmbiguousNucleicAcidIterator(seq)

Base.start(it::AmbiguousNucleicAcidIterator) = find_next_ambiguous(it.seq, 1)
Base.done(it::AmbiguousNucleicAcidIterator, nextpos) = nextpos == 0
function Base.next(it::AmbiguousNucleicAcidIterator, nextpos)
    return nextpos, find_next_ambiguous(it.seq, nextpos + 1)
end

Base.iteratorsize(::AmbiguousNucleicAcidIterator) = Base.SizeUnknown()

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
4x3 Array{Bio.Seq.DNA,2}:
 DNA_A  DNA_A  DNA_A
 DNA_T  DNA_T  DNA_T
 DNA_C  DNA_C  DNA_C
 DNA_G  DNA_G  DNA_G

 julia> seqmatrix(seqs, :seq)
 3x4 Array{Bio.Seq.DNA,2}:
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
    majorityvote{A<:NucleicAcidAlphabet}(seqs::AbstractVector{BioSequence{A}})

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
function majorityvote{A<:NucleicAcidAlphabet}(seqs::AbstractVector{BioSequence{A}})
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
