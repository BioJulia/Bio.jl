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
            nodes[root] |= 1 << 1
            nodes[root] |= idx << 2
            return
        end
    end

    # find five consecutive unused nodes {'#', 'A', 'C', 'G', 'T'}
    base = 0
    for l in root+1:endof(nodes)-4
        if (nodes[l] | nodes[l+1] | nodes[l+2] | nodes[l+3] | nodes[l+4]) & 1 == 0
            base = l
            break
        end
    end
    if base == 0
        resize!(nodes, length(nodes) + 5)
        for l in endof(nodes)-4:endof(nodes)
            nodes[l] = 0
        end
        base = endof(nodes) - 5
    end

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
    nodes[root] = (base << 2) | 1
    for m in 1:5
        nodes[base+m] = 1  # set used bit
    end
    for m in 1:5
        if m == 5
            r = starts[m]:last(range)
        else
            r = starts[m]:starts[m+1]-1
        end
        recursive!(nodes, base + m, seqs, r, j + 1)
    end
end

function Base.in(key, trie::StaticTrie)
    s = 1
    j = 1
    node = trie.nodes[s]
    base = node >> 2
    while node & 3 == 1 && base != 0 && j ≤ endof(key) + 1
        if j == endof(key) + 1
            s = base + 1
        else
            s = base + Int(key[j]) + 2
        end
        node = trie.nodes[s]
        base = node >> 2
        j += 1
    end
    return j == endof(key) + 2 && node & 2 != 0
end
