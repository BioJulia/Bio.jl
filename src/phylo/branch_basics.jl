# phylo/branch_basics.jl
# ======================
#
# Types and methods for phylogenetic trees.
#
# Part of the Bio.Phylo module.
#
# This file contains methods for querying and setting attributes on branches,
# as well as adding and removing branches.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

typealias BranchFields Union{Symbol, AbstractString, Integer}

immutable Branch{T <: BranchFields}
    from::T
    to::T
end


# Functions responsible for getting and setting branch data.
function branchdata{C,B}(tree::Phylogeny{C,B}, edge::Edge)
    return get(tree.edgedata, edge, empty_branch_data(B))
end


function branchdata!{C,B}(tree::Phylogeny{C,B}, edge::Edge, data::B)
    if has_edge(tree.graph, edge)
        tree.edgedata[edge] = data
    else
        warn("Tried to set branch data to a non-existant branch.")
    end
    return tree
end

# Creating and destroying branches in a phylogeny between two nodes.
function add_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge, branchdata::B)
    add_edge!(tree.graph, branch)
    branchdata!(tree, branch, branchdata)
    return tree
end

function add_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge)
    return add_branch!(tree, branch, empty_branch_data(B))
end


function rem_branch!(tree::Phylogeny, branch::Edge)
    rem_edge!(tree.graph, src(branch), dst(branch))
    delete!(tree.edgedata, branch)
    return tree
end

function parent_branch(tree::Phylogeny, vertex::Int)
    @assert hasparent(tree, vertex) error("Vertex is not connected to a parent.")
    return in_edges(tree.graph, vertex)[1]
end

function child_branches(tree::Phylogeny, vertex::Int)
    return out_edges(tree.graph, vertex)
end


# Get and set branchlength mechanism for phylogeneies.

# The basic function on which getting a branch length depends,
# For any kind of phylogeny.
function branchlength(tree::Phylogeny, edge::Edge)
    return branchlength(branchdata(tree, edge))
end

# The basic function on which setting a branch length depends,
# for any kind of phylogeny. The assumption is metadata values on branches are
# immutables or value types - hence the reassignment.
function branchlength!{C,B}(tree::Phylogeny{C,B}, edge::Edge, value::BranchLength)
    branchdata!(tree, edge, branchlength!(branchdata(tree, edge), value))
    return tree
end

function empty_branch_data{T}(::Type{T})
    return T()
end

# To allow branchlength! to work with your metadata type, the type needs the
# following methods defined:
#
# branchlength - which gets the branchlength from the metadata type.
# branchlength!
# optional emptyBranchData - by default it is a no-arg constructor.

empty_branch_data{T<:AbstractFloat}(::Type{T}) = convert(T, -1.0)
branchlength{T<:AbstractFloat}(metadata::T) = metadata
branchlength!{A<:AbstractFloat,B<:AbstractFloat}(metadata::A, bl::B) = convert(A, bl)
