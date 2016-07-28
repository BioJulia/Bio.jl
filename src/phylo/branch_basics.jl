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


"""
    empty_branch_data{T}(::Type{T})

Create a default (empty) instance of any type of branch metadata.

By default this just calls a zero-argument constructor of a type.
But this may be overloaded in some cases.
"""
function empty_branch_data{T}(::Type{T})
    return T()
end

"""
    branchdata{C,B}(tree::Phylogeny{C,B}, edge::Edge)

Getter for metadata associated with branch represented by `edge`.
"""
# Functions responsible for getting and setting branch data.
function branchdata{C,B}(tree::Phylogeny{C,B}, edge::Edge)
    return get(tree.edgedata, edge, empty_branch_data(B))
end

"""
    branchdata!{C,B}(tree::Phylogeny{C,B}, edge::Edge, data::B)

Setter for metadata associated with branch represented by `edge`.
"""
function branchdata!{C,B}(tree::Phylogeny{C,B}, edge::Edge, data::B)
    if has_edge(tree.graph, edge)
        tree.edgedata[edge] = data
    else
        warn("You tried to set branch data to a non-existant branch. Consequently, nothing was done.")
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
    return add_branch!(tree, branch, empty_branch_data(B))
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
function branchlength!{C,B,T<:AbstractFloat}(tree::Phylogeny{C,B}, edge::Edge, value::T)
    branchdata!(tree, edge, branchlength!(branchdata(tree, edge), value))
    return tree
end



# To allow branchlength! to work with your metadata type, the type needs the
# following methods defined:
#
# branchlength - which gets the branchlength from the metadata type.
# branchlength!
# optional empty_branch_data - by default this is a no-arg constructor.

empty_branch_data{T<:AbstractFloat}(::Type{T}) = convert(T, -1.0)
branchlength{T<:AbstractFloat}(metadata::T) = metadata
branchlength!{A<:AbstractFloat,B<:AbstractFloat}(metadata::A, bl::B) = convert(A, bl)
