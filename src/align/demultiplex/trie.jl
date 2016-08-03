# Static Trie
# ===========
#
# Double-array trie for DNA sequences (internal use only).
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable StaticTrie
    nodes::Vector{Int64}
end

function StaticTrie(seqs, values=1:endof(seqs))
    @assert length(seqs) == length(values)
    if !issorted(seqs)
        error("sequences must be sorted")
    end
    for i in 1:endof(seqs)-1
        if seqs[i] == seqs[i+1]
            error("duplicated sequences")
        end
    end
    nodes = Int64[1]
    recursive!(nodes, 1, seqs, values, 1:endof(seqs), 1)
    return StaticTrie(nodes)
end

function recursive!(nodes, root, seqs, values, range, j)
    n_seqs = length(range)
    if n_seqs == 0
        return
    elseif n_seqs == 1
        idx = first(range)
        if j > length(seqs[idx]) + 1
            # leaf
            #nodes[root] = idx << 1
            nodes[root] = values[idx] << 1
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
    while i ≤ last(range) && length(seqs[i]) < j
        i += 1
    end
    for m in 1:4
        starts[m+1] = i
        while i ≤ last(range) && Int(seqs[i][j]) + 1 == m
            i += 1
        end
    end

    # set nodes
    nodes[root] = (base << 1) | 1
    for m in 1:5
        if m == 5
            r = starts[m]:last(range)
        else
            r = starts[m]:starts[m+1]-1
        end
        recursive!(nodes, base + m, seqs, values, r, j + 1)
    end
end

function Base.in(seq, trie::StaticTrie)
    return longest_prefix_length(trie, seq) == length(seq)
end

function findbarcode(trie::StaticTrie, seq)
    s = 1
    for nt in seq
        t = (trie.nodes[s] >> 1) + Int(nt) + 2
        if isleaf(trie.nodes[t])
            break
        end
        s = t
    end
    s = (trie.nodes[s] >> 1) + 1
    return trie.nodes[s] >> 1
end

function isleaf(node::Int64)
    return node & 1 == 0
end

type Demultiplexer
    trie::StaticTrie
    barcodes::Vector{DNASequence}
end

function Demultiplexer(barcodes::Vector{DNASequence}, n_max_mutations::Integer)
    barcodes′, ids = extend(barcodes, n_max_mutations)
    ord = sortperm(barcodes′)
    trie = StaticTrie(barcodes′[ord], ids[ord])
    return Demultiplexer(trie, barcodes)
end

using Iterators

function extend(barcodes, n_max_mutations)
    ret = copy(barcodes)
    ids = collect(1:endof(barcodes))
    for (i, barcode) in enumerate(barcodes)
        if length(barcode) < n_max_mutations
            error("too short barcode")
        end
        for m in 1:n_max_mutations
            circle = hamming_circle(barcode, m)
            append!(ret, circle)
            append!(ids, collect(repeated(i, length(circle))))
        end
    end
    return ret, ids
end

function hamming_circle(seq, m)
    @assert m > 0
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

function demultiplex(x::Demultiplexer, seq::DNASequence)
    return findbarcode(x.trie, seq)
end

#=
function longest_prefix_length(trie::StaticTrie, seq)
    j, _ = traverse(trie, seq)
    return j - 2
end

function traverse(trie::StaticTrie, seq)
    j = 1
    s = 1
    node = trie.nodes[s]
    while !isleaf(node) && j ≤ endof(seq) + 1
        base = node >> 1
        if j == endof(seq) + 1
            s = base + 1
        else
            s = base + Int(seq[j]) + 2
        end
        node = trie.nodes[s]
        j += 1
    end
    return j, s
end
=#
