import Base.@deprecate, Base.@deprecate_binding
import Base.depwarn

# v0.2
# ----

"""
    DNAKmerCounts{k}(seq::Sequence[, step=1])
    RNAKmerCounts{k}(seq::Sequence[, step=1])

Count ocurrences of short (≤ 32) `k`-mers in a sequence.

This method uses a dense table to store counts, requiring O(4^k) memory, so is
not recommended for large k-mer sizes.

# Arguments
* `seq`: A nucleotide sequence.
* `step=1`: K-mers counted are separated by this many nucleotides.
"""
immutable KmerCounts{T,K}
    data::Vector{UInt32}

    function KmerCounts(seq::Sequence, step::Integer=1)
        data = zeros(UInt32, 4^K)
        @inbounds for (_, x) in each(Kmer{T,K}, seq, step)
            data[convert(UInt64, x) + 1] += 1
        end
        return new(data)
    end
end

typealias DNAKmerCounts{K} KmerCounts{DNANucleotide,K}
typealias RNAKmerCounts{K} KmerCounts{RNANucleotide,K}

function Base.getindex{T,K}(counts::KmerCounts{T,K}, x::Kmer{T,K})
    @inbounds c = counts.data[convert(UInt64, x) + 1]
    return c
end

function Base.show{T,K}(io::IO, counts::KmerCounts{T,K})
    println(io, (T == DNANucleotide ? "DNA" : "RNA"), "KmerCounts{", K, "}:")
    for x in UInt64(1):UInt64(4^K)
        s = string(Kmer{T,K}(x - 1))
        println(io, "  ", s, " => ", counts.data[x])
    end
end

@deprecate KmerCounts composition
