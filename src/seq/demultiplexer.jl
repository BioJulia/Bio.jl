# Demultiplexer
# =============
#
# Sequence demultiplexer based on DNA barcodes.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Trie-like data structure for DNA barcodes.
immutable BarcodeTrie
    nodes::Vector{Int64}
end

function BarcodeTrie(barcodes, ids)
    @assert length(barcodes) == length(ids)
    if !issorted(barcodes)
        error("DNA barcodes must be sorted")
    end
    for i in 1:endof(barcodes)-1
        if barcodes[i] == barcodes[i+1]
            error("duplicated sequences")
        end
    end
    nodes = Int64[1]
    build_trie!(nodes, 1, barcodes, ids, 1:endof(barcodes), 1)
    return BarcodeTrie(nodes)
end

# Build a trie recursively.
function build_trie!(nodes, root, barcodes, ids, range, j)
    n = length(range)
    if n == 0
        return
    elseif n == 1
        i = first(range)
        if j > length(barcodes[i]) + 1
            # leaf
            nodes[root] = ids[i] << 1
            return
        end
    end

    # allocate five consecutive nodes {'#', 'A', 'C', 'G', 'T'}
    resize!(nodes, length(nodes) + 5)
    for l in endof(nodes)-4:endof(nodes)
        nodes[l] = 0
    end
    base = endof(nodes) - 5

    # calculate ranges for the next level
    starts = Vector{Int}(5)
    i = first(range)
    starts[1] = i
    while i ≤ last(range) && length(barcodes[i]) < j
        i += 1
    end
    for k in 1:4
        starts[k+1] = i
        while i ≤ last(range) && Int(barcodes[i][j]) + 1 == k
            i += 1
        end
    end

    # set nodes
    nodes[root] = (base << 1) | 1
    for k in 1:5
        if k == 5
            r = starts[k]:last(range)
        else
            r = starts[k]:starts[k+1]-1
        end
        build_trie!(nodes, base + k, barcodes, ids, r, j + 1)
    end
end

# Find the longest match of `seq` prefixes in `trie`.
function findbarcode(trie::BarcodeTrie, seq)
    s = 1
    for nt in seq
        t = (trie.nodes[s] >> 1) + min(Int(nt), 3) + 2
        if trie.nodes[t] & 1 == 0
            # leaf
            break
        end
        s = t
    end
    s = (trie.nodes[s] >> 1) + 1
    return trie.nodes[s] >> 1
end

type Demultiplexer{S}
    # barcode tries with ≥ 0 errors in ascending order
    # (i.e. tries[1]: no errors, tries[2]: one error, ...)
    tries::Vector{BarcodeTrie}

    # original barcodes
    barcodes::Vector{S}

    # distance (:hamming or :levenshtein)
    distance::Symbol
end

function Base.show(io::IO, demultiplexer::Demultiplexer)
    println(io, summary(demultiplexer), ":")
    println(io, "  distance: ", demultiplexer.distance)
    println(io, "  number of barcodes:", length(demultiplexer.barcodes))
      print(io, "  number of correctable errors: ", length(demultiplexer.tries) - 1)
end

"""
    Demultiplexer(barcodes::Vector{DNASequence};
                  n_max_errors::Integer=1,
                  distance::Symbol=:hamming)

Create a demultiplexer object from `barcodes`.

# Arguments
* `barcodes`: the sorted list of DNA barcodes; integer IDs are assigned in this order.
* `n_max_errors=1`: the number of maximum correctable errors.
* `distance=:hamming`: the distance metric (`:hamming` or `:levenshtein`).
"""
function Demultiplexer(barcodes::Vector{DNASequence};
                       n_max_errors::Integer=1,
                       distance::Symbol=:hamming)
    if !issorted(barcodes)
        error("barcodes must be sorted")
    elseif n_max_errors < 0
        error("n_max_errors must be non-negative")
    elseif distance ∉ (:hamming, :levenshtein)
        error("distance must be either :hamming or :levenshtein")
    end

    # check barcodes
    for barcode in barcodes
        for i in 1:endof(barcode)
            if isambiguous(barcode[i])
                error("barcode must be umambiguous")
            end
        end
    end

    # build trie trees for each error level (0-n_max_errors)
    tries = BarcodeTrie[]
    for m in 0:n_max_errors
        # generate "erroneous" barcodes
        barcodes′ = DNASequence[]
        ids = Int[]
        for (i, barcode) in enumerate(barcodes)
            if distance == :hamming
                circle = hamming_circle(barcode, m)
            elseif distance == :levenshtein
                circle = levenshtein_circle(barcode, m)
            else
                # unreachable: already checked above
                assert(false)
            end
            append!(barcodes′, circle)
            append!(ids, collect(repeated(i, length(circle))))
        end
        ord = sortperm(barcodes′)
        push!(tries, BarcodeTrie(barcodes′[ord], ids[ord]))
    end

    return Demultiplexer(tries, barcodes, distance)
end

"""
    demultiplex(demultiplexer::Demultiplexer, seq::Sequence, linear_search_fallback::Bool=false)

Return a barcode index that matches `seq` with least errors.

The upper limit of the number of maximum errors is bounded by the `n_max_errors`
parameter of `demultiplexer`. When `linear_search_fallback` is `true`, this
function tries to find the best matching barcode using linear search.
"""
function demultiplex(demultiplexer::Demultiplexer, seq::Sequence, linear_search_fallback::Bool=false)
    if eltype(seq) != DNANucleotide
        error("sequence must be a DNA sequence")
    end
    for trie in demultiplexer.tries
        i = findbarcode(trie, seq)
        if i != 0
            return i
        end
    end
    if !linear_search_fallback
        # not found
        return 0
    end

    # find the best matching barcode using linear search
    i_min = 0
    dist_min = typemax(Int)
    for (i, barcode) in enumerate(demultiplexer.barcodes)
        if demultiplexer.distance == :hamming
            dist = hamming_distance(barcode, seq[1:length(barcode)])
        elseif demultiplexer.distance == :levenshtein
            dist = sequencelevenshtein_distance(barcode, seq)
        else
            assert(false)
        end
        if dist < dist_min
            dist_min = dist
            i_min = i
        end
    end
    @assert i_min != 0

    return i_min
end

function Base.getindex(demultiplexer::Demultiplexer, i::Integer)
    return demultiplexer.barcodes[i]
end

# Generate a list of sequences s.t. `hamming_distance(seq, seq′) == m`.
function hamming_circle(seq, m)
    if m == 0
        return [seq]
    end
    ret = DNASequence[]
    for ps in combinations(1:endof(seq), m)
        for rs in product(repeated(1:3, m)...)
            seq′ = copy(seq)
            for (p, r) in zip(ps, rs)
                if Int(seq[p]) + 1 ≤ r
                    r += 1
                end
                seq′[p] = (DNA_A:DNA_T)[r]
            end
            push!(ret, seq′)
        end
    end
    return ret
end

# Calculate the Hamming distance between `seq1` and `seq2`.
function hamming_distance(seq1, seq2)
    @assert length(seq1) == length(seq2)
    n = 0
    for (x, y) in zip(seq1, seq2)
        n += x != y
    end
    return n
end

# Generate a list of sequences s.t. `levenshtein_distance(seq, seq′) == m`.
function levenshtein_circle(seq, m)
    if m == 0
        return [seq]
    end
    @assert m > 0
    seqs = DNASequence[]
    # substitution
    append!(seqs, hamming_circle(seq, 1))
    # deletion
    for i in 1:endof(seq)
        push!(seqs, deleteat!(copy(seq), i))
    end
    # insertion
    for i in 1:endof(seq), nt in DNA_A:DNA_T
        if nt != seq[i]
            push!(seqs, insert!(copy(seq), i, nt))
        end
    end
    ball = copy(seqs)
    for seq′ in unique(seqs)
        append!(ball, levenshtein_circle(seq′, m - 1))
    end
    ball = sort!(unique(ball))
    ret = DNASequence[]
    for seq′ in ball
        if levenshtein_distance(seq, seq′) == m
            push!(ret, seq′)
        end
    end
    return ret
end

# Calculate the Levenshtein distance between `seq1` and `seq2`.
function levenshtein_distance(seq1, seq2)
    m = length(seq1)
    n = length(seq2)
    dist = Matrix{Int}(m + 1, n + 1)

    dist[1,1] = 0
    for i in 1:m
        dist[i+1,1] = i
    end

    for j in 1:n
        dist[1,j+1] = j
        for i in 1:m
            dist[i+1,j+1] = min(
                dist[i,j+1] + 1,
                dist[i+1,j] + 1,
                dist[i,j]   + ifelse(seq1[i] == seq2[j], 0, 1))
        end
    end

    return dist[m+1,n+1]
end

# Calculate the Sequence-Levenshtein distance between `seq1` and `seq2`.
# Buschmann, Tilo, and Leonid V. Bystrykh. "Levenshtein error-correcting
# barcodes for multiplexed DNA sequencing." BMC bioinformatics 14.1 (2013): 272.
function sequencelevenshtein_distance(seq1, seq2)
    m = length(seq1)
    n = length(seq2)
    dist = Matrix{Int}(m + 1, n + 1)

    dist[1,1] = 0
    for i in 1:m
        dist[i+1,1] = i
    end

    for j in 1:n
        dist[1,j+1] = j
        for i in 1:m
            dist[i+1,j+1] = min(
                dist[i,j+1] + 1,
                dist[i+1,j] + 1,
                dist[i,j]   + ifelse(seq1[i] == seq2[j], 0, 1))
        end
    end

    mindist = dist[end,end]
    for i in 0:m
        mindist = min(mindist, dist[i+1,end])
    end
    for j in 0:n
        mindist = min(mindist, dist[end,j+1])
    end

    return mindist
end
