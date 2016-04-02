#========================================================#
# Types and methods for the annotation of nodes of trees #
#========================================================#

# Ben J. Ward 2015.
"""
A type-alias: Dictionaries with keys of type `PhyNode` and values of type `T`
are aliased as NodeAnnotations{T}.
"""
typealias NodeAnnotations{T} Dict{PhyNode, T}

"""
TreeAnnotations: A type to allow the assignment of data to clades of a phylogeny.
"""
type TreeAnnotations{T}
    phylogeny::Phylogeny
    annotations::NodeAnnotations{T}
end

"""
Construct a TreeAnnotations object which contains no annotations yet.
"""
function TreeAnnotations{T}(x::Phylogeny, ::Type{T})
    return TreeAnnotations(x, NodeAnnotations{T}())
end

"""
Get the annotation associated with a clade / node by a TreeAnnotatons object.
"""
function Base.getindex{T}(x::TreeAnnotations{T}, index::PhyNode)
    return x.annotations[index]
end

"""
Add or set an annotation to be associated with a clade / node by the TreeAnnotations object.
"""
function Base.setindex!{T}(x::TreeAnnotations{T}, item::T, clade::PhyNode)
    # Should add code here to make sure that clade elected is indeed part of the tree.
    x.annotations[clade] = item
end
