# MinHash
# ========
#
# Functions to generate MinHash sketches of biological sequences
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type MinHashSketch
    sketch::Vector{UInt64}
    kmersize::Int

    function MinHashSketch(sketch::Vector, kmersize::Int)
        length(sketch) > 0 || error("Sketch cannot be empty")
        kmersize > 0 || error("Kmersize must be greater than 0")
        new(sketch, kmersize)
    end
end

Base.getindex(s::MASHSketch, part) = getindex(s.sketch, part)

Base.length(s::MinHashSketch) = length(s.sketch)

Base.size(s::MinHashSketch) = (length(s), s.kmersize)

Base.start(s::MinHashSketch) = start(s.sketch)

Base.next(s::MinHashSketch, state) = next(s.sketch, state)

Base.done(s::MinHashSketch, state) = done(s.sketch, state)

function Base.:(==)(a::MASHSketch, b::MASHSketch)
    a.kmersize == b.kmersize && a.sketch == b.sketch
end


function revcomphash(kmer::Kmer)
    k = minimum((kmer, reverse_complement(kmer)))
    return hash(k)
end


function kmerminhash(seq::BioSequence, kmerset, kmerhashes::Vector{UInt64}, k::Int, s::Int)
    typeof(kmerset) <: Set || typeof(kmerset) <: SortedSet || error("Kmerset must be a `Set` or `SortedSet`")
    length(kmerhashes) <= s || error("Kmerhashes cannot be larger than the set size")

    for kmer in each(DNAKmer{k}, seq)
        if length(kmerhashes) == 0
            if length(kmerset) < s
                push!(kmerset, revcomphash(kmer[2]))
            elseif length(kmerset) == s
                kmerset = SortedSet(kmerset)
                for i in kmerset
                    push!(kmerhashes, pop!(kmerset))
                end
            end
        end

        if length(kmerhashes) == s
            h = revcomphash(kmer[2])
            if h < kmerhashes[end]
                i = searchsortedlast(kmerhashes, h)
                if i == 0 && i != kmerhashes[1]
                    pop!(kmerhashes)
                    unshift!(kmerhashes, h)
                elseif h != kmerhashes[i]
                    pop!(kmerhashes)
                    insert!(kmerhashes, i+1, h)
                end
            end
        end

    end
    return (kmerset, kmerhashes)
end


function minhash(seq::BioSequence, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()
    kmerset, kmerhashes = kmerminhash(seq, kmerset, kmerhashes, k, s)
    return MinHashSketch(kmerhashes, k)
end

function minhash{T<:BioSequence}(seqs::Vector{T}, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()

    for seq in seqs
        kmerset, kmerhashes = kmerminhash(seq, kmerset, kmerhashes, k, s)
    end
    return MinHashSketch(kmerhashes, k)
end

function minhash{T<:BioSequence}(seqs::FASTAReader{T}, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()
    for seq in seqs
        kmerset, kmerhashes = kmerminhash(seq.seq, kmerset, kmerhashes, k, s)
    end
    return MinHashSketch(kmerhashes, k)
end
