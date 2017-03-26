# MinHash
# ========
#
# Functions to generate MinHash sketches of biological sequences
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md
using DataStructures: SortedSet

"""
MinHash Sketch type

MinHash sketches are a sorted set of kmer hashes, containing the smallest `s`
hashes for kmers of a given length. The type contains two parameters:

* .sketch: a sorted set of hashes
* .kmersize: the length of kmers used to generate the sketch
"""
immutable MinHashSketch
    sketch::Vector{UInt64}
    kmersize::Int

    function MinHashSketch(sketch::Vector, kmersize::Integer)
        length(sketch) > 0 || error("Sketch cannot be empty")
        kmersize > 0 || error("Kmersize must be greater than 0")
        return new(sketch, kmersize)
    end
end

Base.getindex(s::MinHashSketch, part) = getindex(s.sketch, part)
Base.length(s::MinHashSketch) = length(s.sketch)
Base.start(s::MinHashSketch) = start(s.sketch)
Base.next(s::MinHashSketch, state) = next(s.sketch, state)
Base.done(s::MinHashSketch, state) = done(s.sketch, state)

function Base.:(==)(a::MinHashSketch, b::MinHashSketch)
    return a.kmersize == b.kmersize && a.sketch == b.sketch
end


function kmerminhash!{k}(::Type{DNAKmer{k}}, seq::BioSequence, s::Integer, kmerhashes::Vector{UInt64})
    # generate first `s` kmers
    iter = each(DNAKmer{k}, seq)
    state = start(iter)
    while length(kmerhashes) < s && !done(iter, state)
        (_, kmer), state = next(iter, state)
        h = hash(canonical(kmer)) # hash lexigraphic minimum of kmer and reverse compliment of kmer
        if h ∉ kmerhashes
            push!(kmerhashes, h)
        end
    end

    sort!(kmerhashes)

    # scan `seq` to make a minhash
    while !done(iter, state)
        (_, kmer), state = next(iter, state)
        h = hash(canonical(kmer)) # hash lexigraphic minimum of kmer and reverse compliment of kmer
        if h < kmerhashes[end]
            i = searchsortedlast(kmerhashes, h)
            if i == 0 && h != kmerhashes[1]
                pop!(kmerhashes)
                unshift!(kmerhashes, h)
            elseif h != kmerhashes[i]
                pop!(kmerhashes)
                insert!(kmerhashes, i+1, h)
            end
        end
    end

    return kmerhashes
end

"""
    minhash(seq, k::Integer, s::Integer)

Generate a MinHash sketch of size `s` for kmers of length `k`.
"""
function minhash(seq::BioSequence, k::Integer, s::Integer)
    kmerhashes = kmerminhash!(DNAKmer{k}, seq, s, UInt64[])
    length(kmerhashes) < s && error("failed to generate enough hashes")

    return MinHashSketch(kmerhashes, k)
end

function minhash{T<:BioSequence}(seqs::Vector{T}, k::Integer, s::Integer)
    kmerhashes = UInt64[]
    for seq in seqs
        kmerminhash!(DNAKmer{k}, seq, s, kmerhashes)
    end

    length(kmerhashes) < s && error("failed to generate enough hashes")
    return MinHashSketch(kmerhashes, k)
end

function minhash(seqs::FASTA.Reader, k::Integer, s::Integer)
    kmerhashes = UInt64[]
    for seq in seqs
        kmerminhash!(DNAKmer{k}, sequence(seq), s, kmerhashes)
    end

    length(kmerhashes) < s && error("failed to generate enough hashes")
    return MinHashSketch(kmerhashes, k)
end
