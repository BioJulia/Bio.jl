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


# The basic function on which getting a branch length depends,
# For any kind of phylogeny.
function branchlength(tree::Phylogeny, edge::Edge)
    return branchlength(tree.edgedata[edge])
end

# The basic function on which setting a branch length depends,
# for any kind of phylogeny. The assumption is metadata values on branches are
# immutables or value types - hence the reassignment.
function branchlength!(tree::Phylogeny, edge, value::Float64)
    tree.edgedata[edge] = branchlength!(tree.edgedata[edge], value)
end


# An example of defining branchlength and branchlength! for a piece of metadata.
# In this case the metadata type is a single Float64, so branches are annotated
# with a branch length and nothing else.
function branchlength(metadata::Float64)
    return metadata
end

function branchlength!(metadata::Float64, value::Float64)
    return value
end


function rem_branch!(tree::Phylogeny, branch::Edge)
    rem_edge!(tree.graph, src(branch), dst(branch))
    delete!(tree.edgedata, branch)
    return tree
end

function add_branch!{C,B}(tree::Phylogeny{C,B}, branch::Edge, branchdata::B)
    add_edge!(tree.graph, branch)
    tree.edgedata[branch] = branchdata
    return tree
end

function parent_branch(tree::Phylogeny, vertex::Int)
    @assert hasparent(tree, vertex) error("Vertex is not connected to a parent.")
    return in_edges(tree.graph, vertex)[1]
end

function child_branches(tree::Phylogeny, vertex::Int)
    return out_edges(tree.graph, vertex)
end
