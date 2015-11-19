# EachKmerIterator
# ================

# Iterate through every k-mer in a nucleotide sequence
immutable EachKmerIterator{T, K}
    seq::NucleotideSequence{T}
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
function each{T, K}(t::Type{Kmer{T, K}}, seq::NucleotideSequence{T}, step::Integer=1)
    @assert K >= 0 "K must be ≥ 0 in EachKmer"
    @assert K <= 32 "K must be ≤ 32 in EachKmer"
    @assert step >= 1 "step must be ≥ 1"

    return EachKmerIterator{T, K}(seq, step)
end

function eachkmer{T}(seq::NucleotideSequence{T}, K::Integer, step::Integer=1)
    return each(Kmer{T, Int(K)}, seq, step)
end

@inline function start{T,K}(it::EachKmerIterator{T,K})
    nextn = findnext(it.seq.ns, it.seq.part.start)
    kmer = Nullable{Tuple{Int,Kmer{T,K}}}()
    kmer, i, nextn = nextkmer(Kmer{T,K}, it.seq, 1, it.step, nextn, kmer)
    return kmer, i, nextn
end

@inline function done{T,K}(it::EachKmerIterator{T,K}, state)
    return isnull(state[1])
end

@inline function next{T,K}(it::EachKmerIterator{T,K}, state)
    kmer, i, nextn = state
    return get(kmer), nextkmer(Kmer{T,K}, it.seq, i, it.step, nextn, kmer)
end

@inline function nextkmer{T,K}(::Type{Kmer{T,K}}, seq, i, step, nextn, kmer)
    # find a position from which we can extract
    # at least K unambiguous nucleotides
    offset = i + seq.part.start - 1
    start = offset
    while nextn > 0 && nextn - offset < K
        i = (nextn + 1) - seq.part.start + 1
        # align
        r = rem(i - 1, step)
        if r > 0
            i += step - r
        end
        offset = i + seq.part.start - 1
        nextn = findnext(seq.ns, offset)
    end
    if seq.part.stop - offset + 1 < K
        # no kmer
        return Nullable{Tuple{Int,Kmer{T,K}}}(), i, nextn
    end

    newkmer = extract_kmer(Kmer{T,K}, seq, offset, kmer)

    # prepare for the next iteration
    if nextn > 0 && nextn < offset + step
        nextn = findnext(seq.ns, offset + step)
    end

    return Nullable((i, newkmer)), i + step, nextn
end

@inline function extract_kmer{K,T}(::Type{Kmer{T,K}}, seq, from, kmer)
    if isnull(kmer) || get(kmer)[1] + K - 1 < from
        # the current kmer doesn't overlap the next kmer
        x = UInt64(0)
        for k in 1:K
            nt = getnuc(T, seq.data, from + k - 1)
            x = x << 2 | UInt8(nt)
        end
    else
        tmp = get(kmer)
        x = let pos = tmp[1], kmer = tmp[2]
            # |<-K->|
            # -------
            #     -------
            # ^   ^
            # pos from
            n = from - pos
            from += K - n
            x = UInt64(kmer)
            for k in 1:n
                nt = getnuc(T, seq.data, from + k - 1)
                x = x << 2 | UInt8(nt)
            end
            x
        end
    end
    return Kmer{T,K}(x)
end
