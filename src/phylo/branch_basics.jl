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

function branchlength{C}(tree::Phylogeny{C, Float64}, edge::Edge)
    return tree.edgedata[edge]
end

function branchlength!{C}(tree::Phylogeny{C, Float64}, edge::Edge, value::Float64)
    tree.edgedata[edge] = value
    return tree
end

function rem_branch!(tree::Phylogeny, branch::Edge)
    rem_edge!(tree.graph, src(branch), dst(branch))
    pop!(tree.edgedata, branch)
    return tree
end

function add_branch!(tree::Phylogeny, branch::Edge, branchlength::Float64)
    add_edge!(tree.graph, branch)
    branchlength!(tree, branch, branchlength)
    return tree
end

function parent_branch(tree::Phylogeny, vertex::Int)
    @assert has_parent(tree, vertex) error("Vertex is not connected to a parent.")
    return in_edges(tree.graph, vertex)[1]
end

function child_branches(tree::Phylogeny, vertex::Int)
    return out_edges(tree.graph, vertex)
end
