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


# Functions responsible for getting and setting branch data.
function branchdata{C,B}(tree::Phylogeny{C,B}, edge::Edge)
    return get(tree.edgedata, edge, B())
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
function add_branch!(tree::Phylogeny, branch::Edge)
    add_edge!(tree.graph, branch)
    return tree
end
function add_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge, branchdata::B)
    add_branch!(tree.graph, branch)
    tree.edgedata[branch] = branchdata
    return tree
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

immutable BranchLength
    value::Float64
end

BranchLength(value::Float64) = BranchLength(value)
BranchLength() = BranchLength(-1.0)
BranchLength(value::BranchLength) = value
BranchLength(metadata::BranchLength, value::BranchLength) = value
Float64(metadata::Float64, value::BranchLength) = value.value

# The basic function on which getting a branch length depends,
# For any kind of phylogeny.
function branchlength(tree::Phylogeny, edge::Edge)
    return BranchLength(branchdata(tree, edge))
end

# The basic function on which setting a branch length depends,
# for any kind of phylogeny. The assumption is metadata values on branches are
# immutables or value types - hence the reassignment.
function branchlength!{C,B}(tree::Phylogeny{C,B}, edge::Edge, value::BranchLength)
    branchdata!(tree, edge, B(branchdata(tree, edge), value))
    return tree
end
