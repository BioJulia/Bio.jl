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

function StaticTrie(seqs)
    if !issorted(seqs)
        error("sequences must be sorted")
    end
    for i in 1:endof(seqs)-1
        if seqs[i] == seqs[i+1]
            error("duplicated sequences")
        end
    end
    nodes = Int64[1]
    recursive!(nodes, 1, seqs, 1:endof(seqs), 1)
    return StaticTrie(nodes)
end

function recursive!(nodes, root, seqs, range, j)
    n_seqs = length(range)
    if n_seqs == 0
        return
    elseif n_seqs == 1
        idx = first(range)
        if j > length(seqs[idx]) + 1
            # leaf
            nodes[root] = idx << 1
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
        recursive!(nodes, base + m, seqs, r, j + 1)
    end
end

function longest_prefix_length(trie::StaticTrie, seq)
    i, node = traverse(trie, seq)
    return i - 2
end

function isleaf(node::Int64)
    return node & 1 == 0
end

function traverse(trie::StaticTrie, seq)
    i = 1
    s = 1
    node = trie.nodes[s]
    while !isleaf(node) && i ≤ endof(seq) + 1
        base = node >> 1
        if i == endof(seq) + 1
            s = base + 1
        else
            s = base + Int(seq[i]) + 2
        end
        node = trie.nodes[s]
        i += 1
    end
    return i, node
end
