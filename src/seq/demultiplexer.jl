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

    # allocate six consecutive nodes {'#', 'A', 'C', 'G', 'T', 'N' (or others)}
    resize!(nodes, length(nodes) + 6)
    for l in endof(nodes)-5:endof(nodes)
        nodes[l] = 0
    end
    base = endof(nodes) - 6

    # calculate ranges for the next level
    starts = Vector{Int}(6)
    i = first(range)
    starts[1] = i  # '#'
    while i ≤ last(range) && length(barcodes[i]) < j
        i += 1
    end
    for k in 1:4  # ACGT
        starts[k+1] = i
        while i ≤ last(range) && findfirst(ACGT, barcodes[i][j]) == k
            i += 1
        end
    end
    starts[6] = i  # N

    # set nodes
    nodes[root] = (base << 1) | 1
    for k in 1:6
        if k == 6
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
        base = trie.nodes[s] >> 1
        if nt == DNA_A
            t = base + 2
        elseif nt == DNA_C
            t = base + 3
        elseif nt == DNA_G
            t = base + 4
        elseif nt == DNA_T
            t = base + 5
        else
            t = base + 6
        end
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
    println(io, "  number of barcodes: ", length(demultiplexer.barcodes))
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
    if n_max_errors < 0
        error("n_max_errors must be non-negative")
    elseif distance ∉ (:hamming, :levenshtein)
        error("distance must be either :hamming or :levenshtein")
    end

    # check barcodes
    for barcode in barcodes
        if hasambiguity(barcode)
            error("barcode must be umambiguous")
        end
    end

    # build trie trees for each error level (0-n_max_errors)
    tries = BarcodeTrie[]
    for m in 0:n_max_errors
        # generate "erroneous" barcodes
        mutated_barcodes = DNASequence[]
        ids = Int[]
        for (i, barcode) in enumerate(barcodes)
            if distance == :hamming
                circle = hamming_circle(barcode, m)
            else
                @assert distance == :levenshtein
                circle = levenshtein_circle(barcode, m)
            end
            append!(mutated_barcodes, circle)
            append!(ids, collect(repeated(i, length(circle))))
        end
        ord = sortperm(mutated_barcodes)
        push!(tries, BarcodeTrie(mutated_barcodes[ord], ids[ord]))
    end

    return Demultiplexer(tries, barcodes, distance)
end

"""
    demultiplex(demultiplexer::Demultiplexer,
                seq::Sequence,
                linear_search_fallback::Bool=false) -> (index, distance)

Return a barcode index that matches `seq` with least errors and its distance.

The upper limit of the number of maximum errors is bounded by the `n_max_errors`
parameter of `demultiplexer`. When `linear_search_fallback` is `true`, this
function tries to find the best matching barcodes using linear search and
returns one of them at random.
"""
function demultiplex(demultiplexer::Demultiplexer, seq::Sequence, linear_search_fallback::Bool=false)
    if eltype(seq) != DNA
        error("sequence must be a DNA sequence")
    end
    for (k, trie) in enumerate(demultiplexer.tries)
        i = findbarcode(trie, seq)
        if i != 0
            return i, k - 1
        end
    end
    if !linear_search_fallback
        # not found
        return 0, -1
    end

    # find the best matching barcode using linear search
    i_min = Int[]
    dist_min = typemax(Int)
    for (i, barcode) in enumerate(demultiplexer.barcodes)
        if demultiplexer.distance == :hamming
            dist = mismatches(barcode, seq[1:length(barcode)])
        elseif demultiplexer.distance == :levenshtein
            dist = sequencelevenshtein_distance(barcode, seq)
        else
            assert(false)
        end
        if dist ≤ dist_min
            if dist < dist_min
                empty!(i_min)
            end
            push!(i_min, i)
            dist_min = dist
        end
    end
    @assert !isempty(i_min)

    return rand(i_min), dist_min
end

function Base.getindex(demultiplexer::Demultiplexer, i::Integer)
    return demultiplexer.barcodes[i]
end

# Generate a list of sequences s.t. `mismatches(seq, seq′) == m`.
function hamming_circle(seq, m)
    if m == 0
        return [seq]
    end
    ret = DNASequence[]
    for ps in Combinatorics.combinations(1:endof(seq), m)
        for rs in Iterators.product(repeated(1:4, m)...)
            seq′ = copy(seq)
            for (p, r) in zip(ps, rs)
                if findfirst(ACGT, seq[p]) ≤ r
                    r += 1
                end
                seq′[p] = ACGTN[r]
            end
            push!(ret, seq′)
        end
    end
    return ret
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
    for i in 1:endof(seq), nt in ACGTN
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
