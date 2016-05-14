# Kmer Iterator
# =============
#
# Iterator over all k-mers in a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Iterate through every k-mer in a nucleotide sequence
immutable EachKmerIterator{T,K,S}
    seq::S
    step::Int
end

# Maybe this function should replace the default constructor.
# Is the (unsafe) default constructor used throughout our code?
"""
Initialize an iterator over all k-mers in a sequence.

Any k-mer containing an N will be skipped over.

### Arguments
  * `t`: Kmer type to enumerate.
  * `seq`: A NucleotideSequence
  * `step`: number of positions between iterated k-mers (default: 1)

### Returns
A EachKmerIterator constructed with these parameters

### Examples

    # iterate over codons
    for (pos, codon) in each(DNAKmer{3}, dna"ATCCTANAGNTACT", 3)
        @show pos, codon
    end
"""
function each{T,K,A<:Union{DNAAlphabet,RNAAlphabet}}(t::Type{Kmer{T,K}}, seq::BioSequence{A}, step::Integer=1)
    @assert K ≥ 0 "K must be ≥ 0 in EachKmer"
    @assert K ≤ 32 "K must be ≤ 32 in EachKmer"
    @assert step ≥ 1 "step must be ≥ 1"
    return EachKmerIterator{T,K,BioSequence{A}}(seq, step)
end

eachkmer{A<:DNAAlphabet}(seq::BioSequence{A}, K::Integer, step::Integer=1) = each(DNAKmer{Int(K)}, seq, step)
eachkmer{A<:RNAAlphabet}(seq::BioSequence{A}, K::Integer, step::Integer=1) = each(RNAKmer{Int(K)}, seq, step)

Base.eltype{T,k,S}(::Type{EachKmerIterator{T,k,S}}) = Tuple{Int,Kmer{T,k}}

if VERSION > v"0.5-"
    Base.iteratorsize(::EachKmerIterator) = Base.SizeUnknown()
end

@inline function Base.start{T,K}(it::EachKmerIterator{T,K})
    nextn = find_next_ambiguous(it.seq, it.seq.part.start)
    pair = Nullable{Tuple{Int,Kmer{T,K}}}()
    pair, nextn = nextkmer(Kmer{T,K}, it.seq, 1, it.step, nextn, pair)
    return pair, nextn
end

@inline function Base.done{T,K}(::EachKmerIterator{T,K}, state)
    return isnull(state[1])
end

@inline function Base.next{T,K}(it::EachKmerIterator{T,K}, state)
    pair, nextn = state
    from, kmer = get(pair)
    return (from, kmer), nextkmer(Kmer{T,K}, it.seq, from + it.step, it.step, nextn, pair)
end

@inline function nextkmer{T,K}(::Type{Kmer{T,K}}, seq, i, step, nextn, kmer)
    # find a position from which we can extract at least K unambiguous nucleotides
    from = i
    while nextn > 0 && nextn - from < K
        from = nextn + 1
        # align `from` since it must be a multiple of `step`
        r = rem(from - 1, step)
        if r > 0
            from += step - r
        end
        nextn = find_next_ambiguous(seq, from)
    end
    if endof(seq) - from + 1 < K
        # no available kmer
        return Nullable{Tuple{Int,Kmer{T,K}}}(), nextn
    end

    newkmer = extract_kmer(Kmer{T,K}, seq, from, kmer)

    # update `nextn` for the next iteration if needed
    if nextn > 0 && nextn < from + step
        nextn = find_next_ambiguous(seq, from)
    end

    return Nullable((from, newkmer)), nextn
end

@inline function extract_kmer{K,T}(::Type{Kmer{T,K}}, seq, from, pair)
    if isnull(pair) || get(pair)[1] + K - 1 < from
        # the last kmer doesn't overlap the extracting one
        x = UInt64(0)
        for k in 1:K
            nt = unsafe_getindex(seq, from + k - 1)
            x = x << 2 | UInt8(nt)
        end
    else
        pos, kmer = get(pair)
        # |<-K->|
        # -------
        #     -------
        # ^   ^
        # pos from
        n = from - pos
        from += K - n
        x = UInt64(kmer)
        for k in 1:n
            nt = unsafe_getindex(seq, from + k - 1)
            x = x << 2 | UInt8(nt)
        end
    end
    return Kmer{T,K}(x)
end
