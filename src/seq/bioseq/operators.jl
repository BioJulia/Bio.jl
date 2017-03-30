# Basic Operators
# ---------------

"""
Count how many nucleotides satisfy a condition (i.e. f(seq[i]) -> true).

The first argument should be a function which accepts a nucleotide as its parameter.
"""
function Base.count(f::Function, seq::BioSequence)
    n = 0
    @inbounds for x in seq
        if f(x)
            n += 1
        end
    end
    return n
end

# Site counting
# -------------

include("site_counting/site_counting.jl")


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
    majorityvote{A<:NucleicAcidAlphabets}(seqs::AbstractVector{BioSequence{A}})

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
function majorityvote{A<:NucleicAcidAlphabets}(seqs::AbstractVector{BioSequence{A}})
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
