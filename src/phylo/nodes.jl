#=================================================================#
# Representing and manipulating phylogenies using connected nodes #
#=================================================================#

# Phylogeny Node Type
# --------------------

"""
PhyNode represents a node or clade of a phylogenetic tree.

It is designed to be flexible and parametric, to enable
the creation of new kinds of nodes with different data and
'abilities', as well as to take advantage of multiple dispatch
and generated functions as much as possible.

PhyNodes connect to parent, and children PhyNodes,
and are parametric to contain any other type as data.

**Fields**

- a reference to its `parent` PhyNode.
- reference to one or more `children`.
- A support value of some type S.
- A branch value of some type B

The branch field is any value that characterises the branch
of the phylogeny that connects the node to its parent.
It could be a single floating point value, in which case it would represent a
branch length. It could also be Void, or Nullable, in which case it might mean
the tree is a cladogram. But users and developers can use any type to represent
a branch.

Similarly, the support field is any value that describes how well this node /
clade is supported in the phylogeny.
Typically this is a single numeric value, a bootstrap value.
But, like the branch field, it is parametric and so any type might be used to
represent the support of a node.

"""
type PhyNode{B,S,M}
    name::ASCIIString
    parent::PhyNode{B,S,M}
    children::Vector{PhyNode{B,S,M}}
    branch::B
    support::S
    metadata::M

    # Internal constructor is a bit ugly, but it has to be in order to achieve
    # incomplete initialisation to make a node point to itself.
    function PhyNode{B,S,M}(name::ASCIIString,
                     branch::B,
                     support::S,
                     nodedata::M,
                     children::Vector{PhyNode{B,S,M}} = PhyNode{B,S,M}[],
                     parent = nothing)
        x = new()
        x.name = name
        x.branch = branch
        x.support = support
        if parent != nothing
            graft!(parent, x)
        else
            x.parent = x
        end
        x.children = PhyNode{B,S,M}[]
        for child in children
            graft!(x, child)
        end
        x.metadata = nodedata
        return x
    end
end


# Outer constructors for Phylogenetic nodes
#-------------------------------------------

PhyNode() = PhyNode{Void, Void, Void}("", nothing, nothing, nothing)

PhyNode(name::ASCIIString) = PhyNode{Void, Void, Void}(name, nothing, nothing, nothing)

PhyNode{B,S}(name::ASCIIString, b::B, s::S) = PhyNode{B,S,Void}(name, b, s, nothing)


# Basic methods, accessing and manipulating fields of individual nodes
#----------------------------------------------------------------------

"""
Get the name of a PhyNode.

**Parameters:**

* `x`: The PhyNode to get the name of.
"""
function name(x::PhyNode)
    return x.name
end


"""
Set the name of a PhyNode.

**Parameters:**

* `x`:    The PhyNode to set the name of.

* `name`: The name to give the PhyNode.
"""
function name!(x::PhyNode, name::ASCIIString)
    x.name = name
end


"""
Test whether the node has a valid (non-null / missing) branch.

**Parameters:**

* `x`:  The PhyNode to test.
"""
function has_branchlength(x::PhyNode)
  error("No defined method for determining whether a $(typeof(x))'s branch length is considered missing.")
end

has_branchlength{S,M}(x::PhyNode{Void,S,M}) = false

has_branchlength{B <: AbstractFloat,S,M}(x::PhyNode{B,S,M}) = x.branch > B(0)

has_branchlength{B,S,M}(x::PhyNode{Nullable{B},S,M}) = !isnull(x.branch)

function has_branchlength{B <: AbstractFloat, S, M}(x::PhyNode{Nullable{B}, S, M})
  !isnull(x.support) && x.support > S(0)
end


"""
Get the branch of a PhyNode.

The branch of a PhyNode characterises it's connection
to its parent.

**Parameters:**

* `x`: The PhyNode to get the branch length of.
"""
branchlength{B,S,M}(x::PhyNode{B,S,M}) = x.branch


"""
Get the branch length of a PhyNode.

**Parameters:**

* `x`: The PhyNode to get the branch length of.
* `replace_unknown`: The value to return if the branchlength is missing/invalid.
"""
function branchlength{B,S,M}(x::PhyNode{B,S,M}, replace_unknown::B)
    return x.branch
end

function branchlength{B,S,M}(x::PhyNode{Nullable{B},S,M}, replace_unknown::B)
    return get(x.branch, replace_unknown)
end

function branchlength{B,S,M}(x::PhyNode{Void,S,M}, replace_unknown::B)
    return replace_unknown
end


"""
Set the branch length of a PhyNode.

**Parameters:**

* `x`:  The PhyNode to set the branchlength of.

* `bl`: The branch length to give the PhyNode, as a Float64 value.
"""
branchlength!{B,S,M}(x::PhyNode{B,S,M}, bl::B) = x.branch = bl

function branchlength!{B,S,M}(x::PhyNode{Nullable{B},S,M}, bl::B)
    if bl > B(0)
        x.branch = Nullable{B}(bl)
    else
        x.branch = Nullable{B}()
    end
end

function branchlength!{B,S,M}(x::PhyNode{Nullable{B},S,M}, bl::Void)
    x.branch = Nullable{B}()
end


"""
Test whether the confidence in the node is known.

**Parameters:**

* `x`:  The PhyNode to test.
"""
function has_support(x::PhyNode)
  error("No defined method for determining whether a $(typeof(x))'s support is considered missing.")
end

has_support{B,S <: AbstractFloat,M}(x::PhyNode{B,S,M}) = x.support >= S(0)

has_support{B,M}(x::PhyNode{B,Void,M}) = false

has_support{B,S,M}(x::PhyNode{B,Nullable{S},M}) = isnull(x.support)

function has_support{B,S <: AbstractFloat,M}(x::PhyNode{B,Nullable{S},M})
  return isnull(x.support) && x.support > S(0)
end


"""
Get the support value of the node.

**Parameters:**

* `x`:  The PhyNode to return the confidence of.
* `replace_unknown`: The value to return if the confidence is null.
"""
support{B,S,M}(x::PhyNode{B,S,M}) = x.support

function support{B,S,M}(x::PhyNode{B,Nullable{S},M}, replace_unknown::S)
    return get(x.support, replace_unknown)
end

function support{B,S,M}(x::PhyNode{B,Void,M}, replace_unknown::S)
    return replace_unknown
end


"""
Set the support value of the node.

**Parameters:**

* `x`:    The PhyNode to set the confidence of.

* `support`: The value of the confidence to be set.
"""
support!{B,S,M}(x::PhyNode{B,S,M}, support::S) = x.support = support

function support!{B,S,M}(x::PhyNode{B,Nullable{S},M}, support::S)
    if support > 0.0
        x.support = Nullable{S}(support)
    else
        x.support = Nullable{S}()
    end
end

function support!{B,S,M}(x::PhyNode{B,Nullable{S},M}, support::Void)
    x.support = Nullable{S}()
end


"""
Get the children of a node.

**Parameters:**

* `x`: The PhyNode to get children for.
"""
function children(x::PhyNode)
    return x.children
end


"""
Add a node to the `children` array of another node.

**Warning:** this method is considered unsafe because it does not build the two-way link between parent and child. If you want to add a child to a node, you should use `graft!()`, which does ensure the two-way link is built.

**Parameters:**

* `parent`: The PhyNode to add a child to.

* `child`:  The PhyNode to add as a child.

"""
function addchild_unsafe!(parent::PhyNode, child::PhyNode)
    if haschild(parent, child)
        error("The child node is already a child of the parent.")
    end
    push!(parent.children, child)
end


"""
Remove a node from the `children` array of another node.

**Warning:** this method is considered unsafe because it does not destroy any two-way link between parent and child. If you want to remove a child from a node, you should use `prune!()`, which does ensure the two-way link is destroyed.

**Parameters:**

* `parent`:

* `child`:
"""
function removechild_unsafe!(parent::PhyNode, child::PhyNode)
    filter!(x -> !(x === child), parent.children)
end


"""
Get the parent of a node.

**Parameters:**

* `x`: The PhyNode to get the parent of.
"""
function parent(x::PhyNode)
    if parentisself(x)
        warn("Node does not have a parent. It is self referential.")
    end
    return x.parent
end


"""
Remove the parent of a `PhyNode` (thus setting the parent property to be self-referential).

**Parameters:**

* `x`: The PhyNode to remove the parent of.
"""
function parent_unsafe!(x::PhyNode)
    x.parent = x
end


"""
Set the parent of a node.

**Warning:** this method is considered unsafe because it does not build the two-way link between parent and child. If you want to add a child to a node, you should use `graft!()`, which does ensure the two-way link is built.

**Parameters:**

* `parent`:  The PhyNode to set as parent.

* `child`:   The PhyNode to set the parent of.
"""
function parent_unsafe!(parent::PhyNode, child::PhyNode)
    child.parent = parent
end


# Methods for printing nodes in various ways
#--------------------------------------------

"""
Basic show method for a PhyNode.
"""
function Base.show(io::IO, n::PhyNode)
    if isempty(n.children)
        print(io, n.name)
    else
        print(io, "(")
        print(io, join(map(string, n.children), ","))
        print(io, ")")
    end
end



# Methods that test nodes for basic relationships and properties
#----------------------------------------------------------------

"""
Test whether a node is empty.

A node is considered "empty", when the name is an empty string,
the branch length and confidence is null, there are no children nodes,
and the node does not have a parent.

**Parameters:**

* `x`: The PhyNode to test.
"""
function Base.isempty(x::PhyNode)
  return x.name == "" && !has_branchlength(x) && !has_support(x) && isunlinked(x)
end


"""
Test whether a node is a leaf.

A node is a leaf, when it has a parent node, but no children nodes.

**Parameters:**

* `x`: The PhyNode to test.
"""
function isleaf(x::PhyNode)
    return hasparent(x) && !haschildren(x)
end


"""
Test whether a node has children.

**Parameters:**

* `x`: The PhyNode to test.
"""
function haschildren(x::PhyNode)
    return length(x.children) > 0
end


"""
Test whether a node is the parent of another specific node.

**Parameters:**

* `parent`: The potential parent PhyNode to test.

* `child`:  The potential child PhyNode to test.
"""
function haschild(parent::PhyNode, child::PhyNode)
    return in(child, parent.children)
end


"""
Test whether a node is its own parent.

**Parameters:**

* `x`: The PhyNode to test.
"""
function parentisself(x::PhyNode)
    return x.parent === x
end


"""
Test whether a node has a parent.

**Parameters:**

* `x`: The PhyNode to test.
"""
function hasparent(x::PhyNode)
    return !parentisself(x)
end


"""
Test whether a node is the root node.

**Parameters:**

* `x`: The PhyNode to test.
"""
function isroot(x::PhyNode)
    return parentisself(x) && haschildren(x)
end


"""
Test whether a node is unlinked, i.e. has no children and no parent.

**Parameters:**

* `x`: The PhyNode to test.
"""
function isunlinked(x::PhyNode)
    return parentisself(x) && !haschildren(x)
end


"""
Test whether a node is linked, i.e. has one or more children and/or a parent.
**Parameters:**

* `x`: The PhyNode to test.
"""
function islinked(x::PhyNode)
    return hasparent(x) || haschildren(x)
end


"""
Test whether a node is internal, i.e. has a parent and one or more children.

**Parameters:**

* `x`: The PhyNode to test.
"""
function isinternal(x::PhyNode)
    return hasparent(x) && haschildren(x)
end


"""
Test whether a node is preterminal i.e. It's children are all leaves.

**Parameters:**

* `x`: The PhyNode to test.
"""
function ispreterminal(x::PhyNode)
    if isleaf(x)
        return false
    end
    return all([isleaf(i) for i in x.children])
end


"""
Test whether a node is semi-preterminal i.e. Some of it's children are leaves, but not all are.

**Parameters:**

* `x`: The PhyNode to test.
"""
function issemipreterminal(x::PhyNode)
    areleaves = [isleaf(i) for i in x.children]
    return any(areleaves) && !all(areleaves)
end


"""
Test whether a node is ancesteral to one or more other nodes.

**Parameters:**

* `posanc`: The PhyNode to test.

* `nodes`: An array of `PhyNode`s that the test node must be ancestral to.
"""
function isancestral(posanc::PhyNode, nodes::Vector{PhyNode})
    desc = descendants(posanc)
    return all([in(node, desc) for node in nodes])
end


"""
Test whether two PhyNodes are equal.

**Parameters:**

* `x`: The left PhyNode to compare.

* `y`: The right PhyNode to compare.
"""
function Base.(:(==))(x::PhyNode, y::PhyNode)
    x.name == y.name &
    isequal(x.branch, y.branch) &
    isequal(x.support, y.support)

end


# Manipulation of node relationships and tree structure,
# including safe pruning and regrafting.
#--------------------------------------------------------

"""
Count the number of children of a node.

**Parameters:**

* `x`: The PhyNode to count the children of.
"""
function countchildren(x::PhyNode)
    return length(x.children)
end


"""
Get the siblings of a node. Included in output is the input node.

**Parameters:**

* `x`: The PhyNode to get siblings for.
"""
function siblings(x::PhyNode)
    if hasparent(x)
        return filter(y -> !(y === x), x.parent.children)
    end
end


"""
Get the descendants of a node.

**Parameters:**

* `x`: The PhyNode to get the descendants of.

"""
function descendants(x::PhyNode)
    if haschildren(x)
        return collect(PhyNode, DepthFirst(x))
    else
        return PhyNode[]
    end
end


"""
Get the terminal descendants of a node. i.e. Nodes that are leaves, which have the input node as an ancestor.

**Parameters:**

* `x`: The PhyNode to get ther terminal descendants of.
"""
function terminaldescendants(x::PhyNode)
    return searchall(DepthFirst(x), isleaf)
end


"""
Get the most recent common ancestor of an array of nodes.

**Parameters:**

* `nodes`:  An array of `PhyNode`s to find the most common ancestor of.
"""
function mrca(nodes::Vector{PhyNode})
    paths = [collect(Tip2Root(i)) for i in nodes]
    convergence = intersect(paths...)
    return convergence[1]
end


"""
Graft a node onto another node, create a parent-child relationship between them.

**Parameters:**

* `parent`: The PhyNode to add a child to.

* `child`:  The PhyNode to add as a child.
"""
function graft!(parent::PhyNode, child::PhyNode)
    if hasparent(child)
        error("This node is already attached to a parent.")
    end
    parent_unsafe!(parent, child)
    addchild_unsafe!(parent, child)
end


"""
Graft a node onto another node, create a parent-child relationship between them, and associatiing a branch length with the relationship.

**Parameters:**

* `parent`:       The PhyNode to add a child to.

* `child`:        The PhyNode to add as a child.

* `branchlength`: The branch length between parent and child.

"""
function graft!(parent::PhyNode, child::PhyNode, branchlength::Float64)
    graft!(parent, child)
    branchlength!(child, branchlength)
end


"""
Graft one or more nodes onto another node, create a parent-child relationship between each of the grafted nodes and the node they are grafted onto.

**Parameters:**

* `parent`:   The PhyNode to add a child to.

* `children`: The array of PhyNodes to add as a child.
"""
function graft!(parent::PhyNode, children::Vector{PhyNode})
    for i in children
        graft!(parent, i)
    end
end


"""
Destroy the relationship between a PhyNode `x` and its parent, returning the PhyNode.

This method cleanly removes the PhyNode `x` from its parent's `children` array, and removes the `parent` reference from the PhyNode `x`. All other fields of the `child` are left intact.

**Parameters:**

* `x`: The PhyNode prune from its parent.
"""
function prune!(x::PhyNode)
    if hasparent(x)
        removechild_unsafe!(x.parent, x)
        parent_unsafe!(x)
    else
        warn("Can't prune this node from its parent, no parent exists for the node.")
    end
    return x
end


"""
Prune a PhyNode from its parent and graft it to another parent.

**Parameters:**

* `prune`:    The PhyNode to remove from its parent.

* `graftto`:  The PhyNode to become the new parent of `prune`.
"""
function pruneregraft!(prune::PhyNode, graftto::PhyNode)
    x = prune!(prune)
    graft!(graftto, x)
end


"""
Prune a PhyNode from its parent and graft it to another parent, setting the branch length.

**Parameters:**

* `prune`:        The PhyNode to remove from its parent.

* `graftto`:      The PhyNode to become the new parent of `prune`.

* `branchlength`: The branch length.
"""
function pruneregraft!(prune::PhyNode, graftto::PhyNode, branchlength::Float64)
    x = prune!(prune)
    graft!(graftto, x, branchlength)
end


"""
Delete a node, destroying the relationships between it and its parent, and it and its children. The children of the node become the children of the node's former parent.

Returns the deleted node.

**Parameter:**

* `x`: The PhyNode to delete.
"""
function Base.delete!(x::PhyNode)
    deleted = prune!(x)
    graft!(parent(deleted), children(deleted))
    return deleted
end


"""
Detach a subtree at a given node.

Returns a new Phylogeny with the detached node as root.

**Parameters:**

* `x`:          The PhyNode to detach.

* `name`:       The name of the new Phylogeny.

* `rooted`:     Whether the detached subtree is rooted.

* `rerootable`: Whether the detached subtree is rerootable.
"""
function detach!(x::PhyNode, name::ASCIIString = "", rooted::Bool = true, rerootable::Bool = true)
    detached = prune!(x)
    return Phylogeny(name, detached, rooted, rerootable)
end


"""
A NodeCache is a dictionary with PhyNodes as keys, and Vectors of some type
as the content. This is used to cache the content of a tree to avoid too much
repeated tree traversal.
"""
typealias NodeCache{T} Dict{PhyNode, Vector{T}}


"""
Cache the contents of the a phylogenetic node.
"""
function cachenodes!{T}(x::PhyNode, vf::Function, store::NodeCache{T})
    for child in x.children
        cachenodes!(child, vf, store)
    end
    store[x] = Vector{T}()
    if haschildren(x)
        for child in x.children
            append!(store[x], store[child])
        end
    else
        push!(store[x], vf(x))
    end
    return store
end


"""
Sort the descendents of a phylogenetic tree.

Descendents are sorted in a consistent way, based on their
names and the names of their descendents properties. This is important for comparisons of nodes and
trees.
"""
function sortdescendents!(x::PhyNode)
    cache = NodeCache{ASCIIString}()
    cachenodes!(cache, x, name)
    for node in DepthFirst(x)
        if haschildren(node)
            sort!(node.children, by = (n) -> string(sort!(cache[n])))
        end
    end
end
