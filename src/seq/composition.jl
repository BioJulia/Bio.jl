# Composition
# ===========
#
# Sequence composition counter.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type Composition{T} <: Associative{T,Int}
    counts::Vector{Int}
end

function Composition{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    counts = zeros(Int, 16)
    @inbounds for x in seq
        counts[reinterpret(UInt8, x) + 1] += 1
    end
    return Composition{eltype(A)}(counts)
end

function Composition{N<:Nucleotide,L}(tuple::NTuple{L,N})
    counts = zeros(Int, 16)
    @inbounds for x in tuple
        counts[reinterpret(UInt8, x) + 1] += 1
    end
    return Composition{N}(counts)
end

function Composition(seq::ReferenceSequence)
    counts = zeros(Int, 16)
    @inbounds for x in seq
        counts[reinterpret(UInt8, x) + 1] += 1
    end
    return Composition{DNANucleotide}(counts)
end

function Composition{T,k}(kmer::Kmer{T,k})
    counts = zeros(Int, 16)
    @inbounds begin
        counts[Int(DNA_A)+1] = count_a(kmer)
        counts[Int(DNA_C)+1] = count_c(kmer)
        counts[Int(DNA_G)+1] = count_g(kmer)
        counts[Int(DNA_T)+1] = count_t(kmer)  # U when T == RNANucleotide
    end
    return Composition{T}(counts)
end

function Composition(seq::AminoAcidSequence)
    counts = zeros(Int, length(alphabet(AminoAcid)))
    @inbounds for x in seq
        counts[reinterpret(UInt8, x) + 1] += 1
    end
    return Composition{AminoAcid}(counts)
end

function Composition{T,k}(iter::EachKmerIterator{T,k})
    counts = zeros(Int, 4^k)
    for (_, x) in iter
        counts[convert(UInt64, x) + 1] += 1
    end
    return Composition{Kmer{T,k}}(counts)
end

"""
    composition(seq | kmer_iter)

Calculate composition of biological symbols in `seq` or k-mers in `kmer_iter`.
"""
function composition(iter::Union{Sequence,EachKmerIterator})
    return Composition(iter)
end

function Base.:(==){T}(x::Composition{T}, y::Composition{T})
    return x.counts == y.counts
end

function Base.length(comp::Composition)
    return length(comp.counts)
end

function Base.start(comp::Composition)
    return UInt64(0)
end

function Base.done(comp::Composition, i)
    return i ≥ length(comp.counts)
end

function Base.next(comp::Composition, i)
    key = convert(keytype(comp), i)
    count = comp.counts[i + 1]
    return (key => count), i + 1
end

function Base.getindex{T}(comp::Composition{T}, x)
    i = convert(UInt64, convert(T, x))
    if !(0 ≤ i < endof(comp.counts))
        throw(KeyError(x))
    end
    return comp.counts[i+1]
end

function Base.copy{T}(comp::Composition{T})
    return Composition{T}(copy(comp.counts))
end

function Base.merge{T}(comp::Composition{T}, other::Composition{T})
    return merge!(copy(comp), other)
end

function Base.merge!{T}(comp::Composition{T}, other::Composition{T})
    @assert length(comp.counts) == length(other.counts)
    for i in 1:endof(comp.counts)
        comp.counts[i] += other.counts[i]
    end
    return comp
end

function majorityvote{T}(comp::Composition{T}, compat = true)
    if compat
        return _mvote_compat(comp)
    else
        return _mvote(comp)
    end
end

function _mvote{T}(comp::Composition{T})
    m = max(comp.counts)
    maxes = convert(Vector{UInt8}, findin(comp.counts, m))
    nucs = Vector{T}(length(maxes))
    @inbounds for i in 1:length(nucs)
        nucs[i] = reinterpret(T, maxes[i] - 1)
    end
    return nucs
end

function _mvote_compat{T}(comp::Composition{T})
    cs = comp.counts
    gapvotes = cs[1]
    a = cs[2]; c = cs[3]; g = cs[5];
    t = cs[9]; m = cs[4]; r = cs[6];
    s = cs[7]; v = cs[8]; w = cs[10];
    y = cs[11]; h = cs[12]; k = cs[13];
    d = cs[14]; b = cs[15]; n = cs[16]
    avotes = a + m + r + w + v + h + d + n
    cvotes = c + m + s + y + v + h + b + n
    gvotes = g + r + s + k + v + d + b + n
    tvotes = t + w + y + k + h + d + b + n
    votes = zeros(Int, 16)
    votes[1] = gapvotes
    votes[2] = avotes
    votes[3] = cvotes
    votes[5] = gvotes
    votes[9] = tvotes
    m = maximum(votes)
    maxes = convert(Vector{UInt8}, findin(votes, m))
    if length(maxes) == 1
        return reinterpret(T, maxes[1] - 0x01)
    else
        nuc = 0x00
        for max in maxes
            nuc = nuc | (max - 0x01)
        end
        return reinterpret(T, nuc)
    end
end

function Base.summary{T}(::Composition{T})
    if T == DNANucleotide
        return "DNA Nucleotide Composition"
    elseif T == RNANucleotide
        return "RNA Nucleotide Composition"
    elseif T == AminoAcid
        return "Amino Acid Composition"
    else
        return string(T, " Composition")
    end
end
