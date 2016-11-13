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


# Any type used as branch metadata should have the following methods defined.

# An empty constructor.
# A constructor that takes a BranchLength type as a

abstract BranchMetaData

immutable BasicBranch <: BranchMetaData
    len::Nullable{Float64}
    conf::Nullable{Float64}
end

@inline function BasicBranch()
    return BasicBranch(Nullable{Float64}(), Nullable{Float64}())
end

@inline function BasicBranch{T}(x::BasicBranch, len::BranchLength{T})
    return BasicBranch(convert(Nullable{Float64}, len), x.conf)
end

function branchlength(x::BasicBranch)
    return x.len
end


"""
    branchdata{C,B}(tree::Phylogeny{C,B}, edge::Edge)

Getter for metadata associated with branch represented by `edge`.
"""
function branchdata{C,B}(tree::Phylogeny{C,B}, edge::Edge)
    return get(tree.edgedata, edge, B())
end

"""
    branchdata!{C,B}(tree::Phylogeny{C,B}, edge::Edge, data::B)

Setter for metadata associated with branch represented by `edge`.
"""
function branchdata!{C,B}(tree::Phylogeny{C,B}, edge::Edge, data::B)
    if has_edge(tree.graph, edge)
        tree.edgedata[edge] = data
    else
        warn("You tried to assign branch data to a non-existant branch. Consequently, nothing was done.")
    end
    return tree
end

"""
    create_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge, branchdata::B)

Create a branch between two nodes in a phylogenetic tree.
Metadata `branchdata` will be associated with the branch.
"""
function create_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge, branchdata::B)
    add_edge!(tree.graph, branch)
    branchdata!(tree, branch, branchdata)
    return tree
end

"""
    create_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge)

Create a branch between two nodes in a phylogenetic tree.
"""
function create_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge)
    return add_branch!(tree, branch, B())
end

"""
    destroy_branch!(tree::Phylogeny, branch::Edge)

Destroy a branch between two nodes in a phylogenetic tree.
When the branch is destroyed so will any associated metadata object.
"""
function destroy_branch!(tree::Phylogeny, branch::Edge)
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

# Lightweight immutable type required to wrap values and distinguish
# them from other values for the purposes of dispatch.
immutable BranchLength{T}
    val::T
end

"""
    branchlength(tree::Phylogeny, edge::Edge)

Get the branchlength of branch defined by `edge` in a phylogeny.
"""
function branchlength(tree::Phylogeny, edge::Edge)
    return branchlength(branchdata(tree, edge))
end

"""
    branchlength!{C,B,T<:AbstractFloat}(tree::Phylogeny{C,B}, edge::Edge, value::T)

Set the branchlength of branch defined by `edge` in a phylogeny.

The assumption is metadata values on branches are immutables or value types,
hence the reassignment using the branchdata! method.
"""
function branchlength!{C,B,T}(tree::Phylogeny{C,B}, edge::Edge, value::T)
    branchdata!(tree, edge, B(branchdata(tree, edge), BranchLength(value)))
    return tree
end
