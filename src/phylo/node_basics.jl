# phylo/node_basics.jl
# ====================
#
# Types and methods for phylogenetic trees.
#
# Part of the Bio.Phylo module.
#
# This file contains the methods for updating or querying characteristics of
# nodes/clades in a phylogeny. In the Phylo representation these are
# vertices of a LightGraphs DiGraph, which are numbered higher than the
# number of taxa + 1.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
    isroot(tree::Phylogeny, vertex::Int)

Test whether a vertex is the root.

# Examples:
```julia
isroot(my_tree, 9)
```
"""
function isroot(tree::Phylogeny, vertex::Int)
    return vertex == root(tree)
end

"""
    isroot{T <: SymString}(tree::Phylogeny, vertex::T)

Test whether a vertex is the root.

# Examples:
```julia
isroot(my_tree, :Root)
isroot(my_tree, :A)
```
"""
function isroot(tree::Phylogeny, name::SymString)
    isroot(tree, tree.vertexindex[name])
end

"""
    isclade(tree::Phylogeny, vertex::Int)

Test whether a node is internal clade, i.e. has a parent and one or more children.

# Examples:
```julia
isclade(my_tree, 6)
```
"""
function isclade(tree::Phylogeny, vertex::Int)
    return vertex >= root(tree)
end

"""
    isclade(tree::Phylogeny, vertex::SymString)

Test whether a node is internal, i.e. has a parent and one or more children.

# Examples:
```julia
isclade(my_tree, :Node7)
isclade(my_tree, "Homo Sapiens")
isclade(my_tree, Symbol("Homo Sapiens"))
```
"""
function isclade(tree::Phylogeny, name::SymString)
    return isclade(tree, tree.vertexindex[name])
end

"""
    isleaf(tree::Phylogeny, vertex::Int)

Test whether a vertex is a leaf.

A vertex is a leaf, when it has a parent node, but no children nodes.

# Examples:
```julia
isleaf(my_phylogeny, 1)
```
"""
function isleaf(tree::Phylogeny, vertex::Int)
    return vertex <= tree.ntaxa
end

"""
    isleaf(tree::Phylogeny, name::SymString)

Test whether a vertex is a leaf.

A vertex is a leaf, when it has a parent node, but no children nodes.

# Examples:
```julia
isleaf(my_phylogeny, :A)
isleaf(my_phylogeny, "A")
```
"""
function isleaf(tree::Phylogeny, name::SymString)
    return isleaf(tree, tree.vertexindex[name])
end

"""
    hasparent(tree::Phylogeny, vertex::Int)

Check if a vertex in a phylogeny is linked to a parent.
"""
function hasparent(tree::Phylogeny, vertex::Int)
    return indegree(tree.graph, vertex) >= 1
end

"""
    hasparent(tree::Phylogeny, name::SymString)

Check if the vertex in a phylogeny is linked to a parent.
"""
function hasparent(tree::Phylogeny, name::SymString)
    return hasparent(tree, tree.vertexindex[name])
end

"""
    hasparent(tree::Phylogeny, vertices::AbstractVector{Int})

Check some vertices in a phylogeny and test if each each is linked to a parent.
"""
function hasparent(tree::Phylogeny, vertices::AbstractVector{Int})
    return indegree(tree.graph, vertices) .>= 1
end

"""
    parent(x::PhyNode, vertex::Int)

Get the parent of a node. In a phylogeny it is not possible
for any one node/vertex to have more than one parent.

# Examples:
```julia
parent(my_tree, 2)
```
"""
function parent(tree::Phylogeny, vertex::Int)
    @assert hasparent(tree, vertex) error("Vertex does not have a parent.")
    return in_neighbors(tree.graph, vertex)[1]
end

"""
    hasparent(tree::Phylogeny, child::Int, parent::Int)

Test if a vertex of a phylogeny is the parent of another vertex.
"""
function hasparent(tree::Phylogeny, child::Int, parent::Int)
    return has_edge(tree.graph, parent, child)
end

"""
    haschildren(tree::Phylogeny, vertex::Int)

Test whether a vertex in a phylogeny has children.
"""
function haschildren(tree::Phylogeny, vertex::Int)
    return nchildren(tree, vertex) > 0
end

"""
    haschildren(tree::Phylogeny, vertices::AbstractVector{Int})

Test vertices in a phylogeny to see if it has children.
"""
function haschildren(tree::Phylogeny, vertices::AbstractVector{Int})
    return nchildren(tree, vertices) .> 0
end

"""
    nchildren(tree::Phylogeny, vertex::Int)

Count the number of clades that are direct children of a clade.
"""
function nchildren(tree::Phylogeny, vertex::Int)
    return outdegree(tree.graph, vertex)
end

"""
    nchildren(tree::Phylogeny, vertices::AbstractVector{Int})

Count the number of clades that are direct children for a selection of clades.
"""
function nchildren(tree::Phylogeny, vertices::AbstractVector{Int})
    return outdegree(tree.graph, vertices)
end

"""
    nchildren(tree::Phylogeny, name::T)

Count the number of clades that are direct children of a clade.
"""
function nchildren(tree::Phylogeny, name::SymString)
    return nchildren(tree, tree.vertexindex[name])
end

"""
    haschild(tree::Phylogeny, child::Int, parent::Int)

Test if a vertex of a phylogeny is the child of another vertex.
"""
function haschild(tree::Phylogeny, parent::Int, child::Int)
    return has_edge(tree.graph, parent, child)
end

"""
    children(tree::Phylogeny, vertex::Int)

Get the children vertices of a vertex in the phylogeny.
"""
function children(tree::Phylogeny, vertex::Int)
    return copy(out_neighbors(tree.graph, vertex))
end

"""
    ispreterminal(tree::Phylogeny)

Test whether a node is preterminal i.e. It's children are all leaves.
"""
function ispreterminal(tree::Phylogeny, vertex::Int)
    if !haschildren(tree, vertex)
        return false
    end
    # We use out_neighbours below as it returns a reference, rather than a copy
    # as the children method above does, so we save unnessecery copy ops.
    return all([isleaf(tree, i) for i in out_neighbors(tree.graph, vertex)])
end

"""
    issemipreterminal(tree::Phylogeny, vertex::Int)

Test whether a node is semi-preterminal i.e. Some of it's children are leaves,
but not all are.
"""
function issemipreterminal(tree::Phylogeny, vertex::Int)
    if !haschildren(tree, vertex)
        return false
    end
    # We use out_neighbours below as it returns a reference, rather than a copy
    # as the children method above does, so we save unnessecery copy ops.
    areleaves = [isleaf(tree, i) for i in out_neighbors(tree.graph, vertex)]
    return any(areleaves) && !all(areleaves)
end

"""
    isredundant(tree::Phylogeny, vertex::Int)

Test if a vertex in a phylogeny is redundant or not.
A vertex is considered redundant is it has a parent and only one child.
"""
function isredundant(tree::Phylogeny, vertex::Int)
    return hasparent(tree, vertex) && nchildren(tree, vertex) == 1
end

"""
    isredundant(tree::Phylogeny, vertices::AbstractArray{Int})

Test if vertices in a phylogeny to see if each is redundant or not.
A vertex is considered redundant is it has a parent and only one child.
"""
function isredundant(tree::Phylogeny, vertices::AbstractArray{Int})
    return hasparent(tree, vertices) & nchildren(tree, vertices) .== 1
end
