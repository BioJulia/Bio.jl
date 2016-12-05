# phylo/manipulation.jl
# =====================
#
# Types and methods for phylogenetic trees.
#
# Part of the Bio.Phylo module.
#
# This file contains methods for manipulating the structure of phylogenies.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md



# Internal functions

"""
    unconnected_clades(tree::Phylogeny)

**Unexported function for internal use.**

Gets the internal clade vertices which are not connected to either parents or
children.
"""
function unconnected_clades(tree::Phylogeny)
    cv = clades(tree)
    idx = find(degree(tree.graph, cv) .== 0)
    return cv[idx]
end

"""
    subtree_roots(tree::Phylogeny)

**Unexported function for internal use.**

Gets clade vertices from the phylogeny which are roots of detached subphylogenies.
Such vertices have children, but no parent, and are not THE root of the
phylogeny (n + 1). Such subtrees are often created after pruning of trees.
"""
function subtree_roots(tree::Phylogeny)
    cv = clades(tree)
    idx = find(haschildren(tree, cv) & !hasparent(tree, cv))
    return cv[idx]
end

"""
    disconnect_root(tree::Phylogeny)

**Unexported function for internal use.**

Unconnects the root vertex from all of its children.
"""
function disconnect_root!(tree::Phylogeny)
    r = root(tree)
    for i in children(tree, r)
        destroy_branch!(tree, Edge(r, i))
    end
    return tree
end

"""
    delete!(tree::Phylogeny, vertex::Int, preserve_bl::Bool = false)

Delete a node from a phylogenetic tree.
"""
function Base.delete!(tree::Phylogeny, vertex::Int, preserve_bl::Bool = false)
    p = Phylo.parent(tree, vertex)
    # Delete the connection to parent but remember the branchlength of the
    # deleted branch.
    parentedge = Edge(p, vertex)
    lentoparent = preserve_bl ? branchlength(tree, parentedge) : 0.0
    destroy_branch!(tree, parentedge)
    for branch in child_branches(tree, vertex)
        newbl = branchlength(tree, branch) + lentoparent
        add_branch!(tree, Edge(parent, dst(branch)), newbl)
        destroy_branch!(tree, Edge(vertex, dst(branch)))
    end
    return tree
end
