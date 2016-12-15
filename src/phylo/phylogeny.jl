# phylo/phylogeny.jl
# ==================
#
# Types and methods for phylogenetic trees.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md



# The Phylogeny Type
# -------------------

# The internal structure and representation of the Phylogeny type uses a DiGraph.
# The vertices of the DiGraph are reserved for different specific purposes.

# Vertex ranges:
# * 1:tree.ntaxa => leaf vertices
# * tree.ntaxa + 1 => root vertex
# * tree.ntaxa + 1:<#nodes> => clade vertices

"""
Phylogeny represents a phylogenetic tree.

The type is parametric with two parameters `C` and `B`.

This is because it is common to want to annotate tips, clades, and branches
in a phylogeny with data to create a richer model of evolution of do other things
like dictate aesthetic values when plotting.

Type parameter `C` dictates what datatype can be stored in the phylogeny to annotate
clades and tips. Type parameter `B` dictates what datatype can be stored in the
phylogeny to annotate branches. Think `C` for clades and `B` for branches.
"""
type Phylogeny{C, B}
    graph::DiGraph
    vertexindex::Indexer{Int}
    vertexdata::Vector{C}
    edgedata::Dict{Edge, B}
    ntaxa::Int
    rooted::Bool
    rerootable::Bool
end

# Phylogeny Constructors
# -----------------------

"""
    Phylogeny{C, B}(::Type{C}, ::Type{B}, tips::Vector{Symbol})

Construct a completely unresolved phylogeny, with no annotated data.

The phylogeny will store annotation data of type `C` for tree nodes
(clades and tips). It will also store data of type `B` associated with branches.

# Examples
```julia
Phylogeny(Void, Void, [:Human, :Chimp, :Dog, :Cat])
Phylogeny(Float64, Float64, [:Human, :Chimp, :Dog, :Cat])
```
"""
function Phylogeny{C,B}(::Type{C}, ::Type{B}, taxa::Vector{Symbol})
    ntaxa = length(taxa)
    maxVertices = maxvertices(ntaxa)
    taxaAndRoot = push!(copy(taxa), :Root)
    idxer = Indexer(taxaAndRoot, Int)
    g = DiGraph(maxVertices)
    for i in 1:ntaxa
        add_edge!(g, ntaxa + 1, i)
    end
    v = Vector{C}()
    d = Dict{Edge,B}()
    return Phylogeny{C,B}(g, idxer, v, d, ntaxa, false, true)
end

"""
    Phylogeny(taxa::Vector{Symbol})

Construct a completely unresolved phylogeny, with no annotated data.

This constructor assumes you want a phylogeny where the only data associated
with it is branch lengths (`Float64`) for the edges, and confidence values
(`Float64`) for the clades.

# Examples
```julia
Phylogeny([:Human, :Chimp, :Dog, :Cat])
```
"""
function Phylogeny(taxa::Vector{Symbol})
    return Phylogeny(Float64, BasicBranch, taxa)
end

"""
    Phylogeny{T <: AbstractString}(taxa::Vector{T})
Construct a completely unresolved phylogeny, with no annotated data.

This constructor assumes you want a phylogeny where the only data associated
with it is branch lengths (`Float64`) for the edges, and confidence values
(`Float64`) for the clades.
"""
function Phylogeny{T <: AbstractString}(taxa::Vector{T})
    return Phylogeny(convert(Vector{Symbol}, taxa))
end

function Base.copy{C, B}(tree::Phylogeny{C, B})
    g = copy(tree.graph)
    i = copy(tree.vertexindex)
    v = copy(tree.vertexdata)
    d = copy(tree.edgedata)
    return Phylogeny{C, B}(g, i, v, d, tree.ntaxa, tree.rooted, tree.rerootable)
end

function Base.:(==)(x::Phylogeny, y::Phylogeny)
    g = x.graph == x.graph
    vi = x.vertexindex == y.vertexindex
    vd = x.vertexdata == y.vertexdata
    ed = x.edgedata == y.edgedata
    t = x.ntaxa == y.ntaxa
    r = x.rooted == y.rooted
    rr = x.rerootable == y.rerootable
    return g && vi && vd && ed && t && r && rr
end

# Miscellaneous Functions
# ------------------------

"""
    n_possible_rooted(n::Int)

Returns the number of possible rooted tree topologies for n taxa.

**Warning:** This method currently overflows when `n` > 12.
"""
function n_possible_rooted(n::Int)
    @assert n >= 2 DomainError("n must be greater or equal to 2.")
    num = factorial(2 * n - 3)
    den = (2^(n - 2)) * factorial(n - 2)
    return Int(num / den)
end

"""
    n_possible_unrooted(n::Int)

Returns the number of possible rooted tree topologies for n taxa.

**Warning:** This method currently overflows when `n` > 13.
"""
function n_possible_unrooted(n::Int)
    @assert n >= 3 DomainError("n must be greater of equal to 3.")
    return n_possible_rooted(n - 1)
end

maxvertices(n::Int) = 2 * n - 1

typealias SymString Union{Symbol, AbstractString}

"""
    isrooted(x::Phylogeny)

Test whether a Phylogeny is rooted.

# Examples
```julia
isrooted(my_phylogeny)
```
"""
function isrooted(x::Phylogeny)
    return x.rooted
end

"""
    isrerootable(x::Phylogeny)

Test whether a Phylogeny is re-rootable.

# Examples
```julia
isrerootable(my_phylogeny)
```
"""
function isrerootable(x::Phylogeny)
    return x.rerootable
end

"Get the vertex of the tree which represents the root of the tree."
root(tree::Phylogeny) = tree.ntaxa + 1

"Get the range of vertices which represent the internal clades of the tree."
clades(tree::Phylogeny) = root(tree):nv(tree.graph)

"Get the range of vertices which represent the leaves of the tree."
leaves(tree::Phylogeny) = 1:tree.ntaxa

"""
    hasparent(tree::Phylogeny)

Check every vertex in a phylogeny and test if each is linked to a parent.
"""
function hasparent(tree::Phylogeny)
    return indegree(tree.graph) .>= 1
end

"""
    haschildren(tree::Phylogeny)

Test each vertex in a phylogeny to see if it has children.
"""
function haschildren(tree::Phylogeny)
    return nchildren(tree) .> 0
end

"""
    nchildren(tree::Phylogeny)

Count the number of clades that are direct children of each clade.
"""
function nchildren(tree::Phylogeny)
    return outdegree(tree.graph)
end

"""
    isredundant(tree::Phylogeny)

Test each vertex in a phylogeny and test if it redundant.

Vertices are considered redundant if they have a parent and 1 child.
"""
function isredundant(tree::Phylogeny)
    return hasparent(tree) & nchildren(tree) .== 1
end
