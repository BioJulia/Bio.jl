# Kmer Iterator
# =============
#
# Iterator over all k-mers in a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Iterate through every k-mer in a nucleotide sequence
immutable EachKmerIterator{T<:Kmer,S<:Sequence}
    seq::S
    step::Int
    start::Int
end

"""
    each(::Type{Kmer{T,k}}, seq::Sequence[, step=1])

Initialize an iterator over all k-mers in a sequence `seq` skipping ambiguous
nucleotides without changing the reading frame.

# Arguments
* `Kmer{T,k}`: k-mer type to enumerate.
* `seq`: a nucleotide sequence.
* `step=1`: the number of positions between iterated k-mers

# Examples
```
# iterate over DNA codons
for (pos, codon) in each(DNAKmer{3}, dna"ATCCTANAGNTACT", 3)
    @show pos, codon
end
```
"""
function each{T,K}(::Type{Kmer{T,K}}, seq::Sequence, step::Integer=1)
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    elseif !(0 ≤ K ≤ 32)
        throw(ArgumentError("k-mer length must be between 0 and 32"))
    elseif step < 1
        throw(ArgumentError("step size must be positive"))
    end
    return EachKmerIterator{Kmer{T,K},typeof(seq)}(seq, step, 1)
end

eachkmer{A<:DNAAlphabet}(seq::BioSequence{A}, K::Integer, step::Integer=1) = each(DNAKmer{Int(K)}, seq, step)
eachkmer{A<:RNAAlphabet}(seq::BioSequence{A}, K::Integer, step::Integer=1) = each(RNAKmer{Int(K)}, seq, step)
eachkmer(seq::ReferenceSequence, K::Integer, step::Integer=1) = each(DNAKmer{Int(K)}, seq, step)

Base.eltype{T,S}(::Type{EachKmerIterator{T,S}}) = Tuple{Int,T}

function Base.iteratorsize{T,S}(::Type{EachKmerIterator{T,S}})
    return Base.SizeUnknown()
end

function Base.start{T}(it::EachKmerIterator{T})
    k = kmersize(T)
    pos = it.start
    kmer::UInt64 = 0
    while pos + k - 1 ≤ endof(it.seq)
        kmer, ok = extract_kmer_impl(it.seq, pos, k)
        if ok
            break
        end
        pos += it.step
    end
    return pos, kmer
end

@inline function Base.done{T}(it::EachKmerIterator{T}, state)
    return state[1] + kmersize(T) - 1 > endof(it.seq)
end

@inline function Base.next{T}(it::EachKmerIterator{T}, state)
    k = kmersize(T)
    pos, kmer = state
    pos += it.step
    has_ambiguous = false

    # faster path: return the next overlapping kmer if possible
    if it.step < k && pos + k - 1 ≤ endof(it.seq)
        offset = k - it.step
        if it.step == 1
            nt = inbounds_getindex(it.seq, pos+offset)
            kmer = kmer << 2 | trailing_zeros(nt)
            has_ambiguous |= isambiguous(nt)
        else
            for i in 1:it.step
                nt = inbounds_getindex(it.seq, pos+i-1+offset)
                kmer = kmer << 2 | trailing_zeros(nt)
                has_ambiguous |= isambiguous(nt)
            end
        end
        if !has_ambiguous
            return (state[1], convert(T, state[2])), (pos, kmer)
        end
    end
    
    # fallback
    while pos + k - 1 ≤ endof(it.seq)
        kmer, ok = extract_kmer_impl(it.seq, pos, k)
        if ok
            break
        end
        pos += it.step
    end
    return (state[1], convert(T, state[2])), (pos, kmer)
end

function extract_kmer_impl(seq, from, k)
    kmer::UInt64 = 0
    has_ambiguous = false
    for i in 1:k
        nt = inbounds_getindex(seq, from+i-1)
        kmer = kmer << 2 | trailing_zeros(nt)
        has_ambiguous |= isambiguous(nt)
    end
    return kmer, !has_ambiguous
end
