"""
Phylogeny represents a phylogenetic tree.

A tree can have:

- `name`
- `root`
- `rooted`
- `rerootable`
"""
type Phylogeny
    name::ASCIIString
    root::PhyNode
    rooted::Bool
    rerootable::Bool

    Phylogeny() = new("", PhyNode(), false, true)
end


"""
Create a Phylogeny with a name, root node, and set whether it is rooted and whether
it is re-rootable.

**Parameters:**
* `name`:       The name of the tree.
* `root`:       The root node.
* `rooted`:     Whether the tree is rooted.
* `rerootable`: Whether the tree is re-rootable.
"""
function Phylogeny(name::ASCIIString, root::PhyNode, rooted::Bool, rerootable::Bool)
    x = Phylogeny()
    name!(x, name)
    x.root = root
    x.rooted = rooted
    rerootable!(x, rerootable)
    return x
end


"""
Test whether a phylogeny is empty.

**Parameters:**
* `x`: The Phylogeny to test.
"""
function isempty(x::Phylogeny)
    return isempty(x.root)
end


"""
Set the name of a Phylogeny

**Parameters:**
* `x`:    The Phylogeny to set the name of.
* `name`: The name to set.
"""
function name!(x::Phylogeny, name::ASCIIString)
    x.name = name
end


"""
Test whether a Phylogeny is rooted.

**Parameters:**
* `x`: The Phylogeny to test.
"""
function isrooted(x::Phylogeny)
    return x.rooted
end


"""
Test whether a Phylogeny is re-rootable.

**Parameters:**
* `x`: The Phylogeny to test.
"""
function isrerootable(x::Phylogeny)
    return x.rerootable
end


"""
Get the root node of a Phylogeny.

**Parameters:**
* `x`: The Phylogeny to get the root of.
"""
function root(x::Phylogeny)
    return x.root
end


"""
Test whether a given node is in a given tree.

**Parameters:**
* `tree`:   The Phylogeny to check.
* `clade`:  The PhyNode to check.
"""
function isintree(tree::Phylogeny, clade::PhyNode)
    s = search(BreadthFirst(tree), x -> x === clade)
    return typeof(s) == PhyNode
end


"""
Root a tree at the midpoint between the two most distant taxa.

This method modifies the `tree` variable.

**Parameters:**
* `tree`: The Phylogeny to root.
"""
function root!(tree::Phylogeny, newbl::Float64 = -1.0)
    midpoint = findmidpoint(tree)
    root!(tree, midpoint, newbl)
end


"""
Find the maximum branch length in a dictionary mapping nodes to their branch lengths.

**Parameters:**
* `dict`: The dictionary.
"""
function maxindict(dictionary::Dict)
    keyvalpairs = collect(dictionary)
    values = [i[2] for i in keyvalpairs]
    matches = maximum(values) .== values
    return keyvalpairs[matches][1]
end


"""
Find the node that is furthest from the root of a tree.

**Parameters:**
* `tree` The Phylogeny to search.
"""
function furthestfromroot(tree::Phylogeny)
    distances = distance(tree)
    return maxindict(distances.annotations)
end


"""
Find the leaf that is furthest from a given node in a tree.

  **Parameters:**
  * `tree`: The Phylogeny containing the nodes.
  * `node`: The PhyNode find the furthest node from.
"""
function furthestleaf(tree::Phylogeny, node::PhyNode)
    distances = Dict{PhyNode, Float64}([i => distance(tree, node, i) for i in terminaldescendents(root(tree))])
    return maxindict(distances)
end


"""
Find the midpoint of a tree.

**Parameters:**
* `tree`: The Phylogeny to find the midpoint of.
"""
function findmidpoint(tree::Phylogeny)
    furthestfromroot, ffrdist = furthestfromroot(tree)
    furthestfromleaf, ffldist = furthestleaf(tree, furthestfromroot)
    outgroup = furthestfromroot
    middistance = ffldist / 2.0
    cdist = 0.0
    current = furthestfromroot
    while true
        cdist += branchlength(current)
        if cdist > middistance
            break
        else
            current = parent(current)
        end
    end
    return current
end


"""
Root a tree using a given array of nodes as the outgroup, and optionally setting the branch length.

**Parameters:**
* `tree`: The Phylogeny to root.
* `outgroup`: An array of PhyNodes to use as outgroup.
* `newbl`: The new branch length (optional).
"""
function root!(tree::Phylogeny,
               outgroup::Vector{PhyNode},
               newbl::Float64 = -1.0)
    o = mrca(outgroup)
    root!(tree, o, newbl)
end


"""
Root a tree using a given node as the outgroup, and optionally setting the branch length,

**Parameters:**
* `tree`:     The Phylogeny to root.
* `outgroup`: A PhyNode to use as outgroup.
* `newbl`:    The new branch length (optional).
"""
function root!(tree::Phylogeny, outgroup::PhyNode, newbl::Float64 = -1.0)
    # Check for errors and edge cases first as much as possible.
    # 1 - The tree is not rerootable.
    if !isrerootable(tree)
        error("Phylogeny is not rerootable!")
    end
    # 2 - The specified outgroup is already the root.
    if isroot(outgroup)
        error("New root is already the root!")
    end
    # 3 - Check the new branch length for the outgroup
    # is between 0.0 and the old previous branchlength.
    previousbranchlength = branchlength(outgroup)
    @assert 0.0 <= newbl <= previousbranchlength
    # 4 - Check that the proposed outgroup is indeed part of the tree.
    if !isintree(tree, outgroup)
        error("The specified outgroup is not part of the phylogeny.")
    end

    #  the path from the outgroup to the root, excluding the root.
    outgrouppath = collect(Tip2Root(outgroup))[2:end - 1]

    # Edge case, the outgroup to be the new root
    # is terminal or the new branch length is not nothing,
    # we need a new root with a branch to the outgroup.
    if isleaf(outgroup) || newbl != 0.0
        newroot = PhyNode("NewRoot", branchlength(root(tree)))
        pruneregraft!(outgroup, newroot, newbl)
        if length(outgrouppath) == 0
            # There aren't any nodes between the outgroup
            # and origional group to rearrange.
            newparent = newroot
        else
            parent = splice!(outgrouppath, 1)
            previousbranchlength, parent.branchlength = parent.branchlength, previousbranchlength - branchlength(outgroup)
            pruneregraft!(parent, newroot)
            newparent = parent
        end
    else
        # Use the provided outgroup as a
        # trifurcating root if the node is not a leaf / newbl is 0.0.
        newroot = newparent = outgroup
        branchlength!(newroot, branchlength(root(tree)))
    end

    # Now we trace the outgroup lineage back,
    # reattaching the subclades under the new root!
    for parent in outgrouppath
        #prune!(newparent)
        previousbranchlength, parent.branchlength =
        parent.branchlength, previousbranchlength
        pruneregraft!(parent, newparent)
        newparent = parent
    end

    # Now we have two s of connected PhyNodes.
    # One begins the with the new root and contains the
    # nodes rearranged as per the backtracking process
    # along outgrouppath. The other is the nodes still
    # connected to the old root.
    # This needs to be resolved.

    # If the old root only has one child, it was bifurcating,
    # and if so, must be removed and the branch lengths resolved,
    # appropriately.
    if countchildren(tree.root) == 1
        ingroup = children(root(tree))[1]
        branchlength!(ingroup, branchlength(ingroup) + previousbranchlength)
        pruneregraft!(ingroup, newparent)
    else
        # If the root has more than one child,
        # then it needs to be kept as an internal node.
        branchlength!(tree.root, previousbranchlength)
        graft!(newparent, tree.root)
    end

    # TODO / FUTURE IMPROVEMENT - COPYING OF OLD ROOT ATTRIBUTES OR DATA TO NEW ROOT.

    tree.root = newroot
    tree.rooted = true
end


# This is probably unnecessary given root puts the rooted flag to true.
# perhaps and unroot! method is more appropriate.
"""
Unroot a tree.

**Parameters:**
* `x`: The Phylogeny to unroot.
"""
function unroot!(x::Phylogeny)
    x.rooted = false
end


"""
Set whether a tree is re-rootable.

**Parameters:**
* `x`:          The Phylogeny.
* `rerootable`: Whether the Phylogeny is re-rootable.
"""
function rerootable!(x::Phylogeny, rerootable::Bool)
    x.rerootable = rerootable
end


"""
Get the terminal nodes of a phylogeny.

**Parameters:**
* `x`: The Phylogeny.
"""
function terminals(x::Phylogeny)
    return terminaldescendents(x.root)
end


"""
Get one or more nodes by name.

**Parameters:**
* `tree`:  The Phylogeny to search.
* `names`: The names of the nodes to get.
"""
function getindex(tree::Phylogeny, names::Array{ASCIIString, 1})
  return searchall(DepthFirst(tree), x -> in(name(x), names))
end


"""
Get one node by name.

**Parameters:**
* `tree`:  The Phylogeny to search.
* `names`: The name of the nodes to get.
"""
function getindex(tree::Phylogeny, name::ASCIIString)
    for i in DepthFirst(tree)
        if i.name == name
            return i
        end
    end
    error("No Node in phylogeny by specified name.")
end


"""
Generate an index mapping names to nodes

**Parameters:**
* `tree`: The Phylogeny to index.
"""
function generateindex(tree::Phylogeny)
    output = Dict{ASCIIString, PhyNode}()
    for i = BreadthFirst(tree)
        if haskey(output, name(i))
            error("You are trying to build an index dict " *
                "of a tree with clades of the same name.")
        end
        output[name(i)] = i
    end
    return output
end


"""
Find the shortest path between two nodes in a tree.

**Parameters:**
* `tree`: The Phylogeny to search in .
* `n1`:   The first node.
* `n2`:   The second node.
"""
function pathbetween(tree::Phylogeny, n1::PhyNode, n2::PhyNode)
    if !isintree(tree, n1) || !isintree(tree, n2)
        error("One of the nodes is not present in the tree.")
    end
    p1::Vector{PhyNode} = collect(Tip2Root(n1))
    p2::Vector{PhyNode} = collect(Tip2Root(n2))
    inter::Vector{PhyNode} = intersect(p1, p2)
    filter!((x) -> !in(x, inter), p1)
    filter!((x) -> !in(x, inter), p2)
    return [p1, inter[1], reverse(p2)]
end


"""
Find the distance between two nodes in a tree.

**Parameters:**
* `tree`: The Phylogeny to search in.
* `n1`:   The first node.
* `n2`:   The second node.
"""
function distance(tree::Phylogeny, n1::PhyNode, n2::PhyNode)
    if !isintree(tree, n1) || !isintree(tree, n2)
        error("One of the nodes is not present in the tree.")
    end
    p1::Vector{PhyNode} = collect(Tip2Root(n1))
    p2::Vector{PhyNode} = collect(Tip2Root(n2))
    inter::Vector{PhyNode} = intersect(p1, p2)
    filter!((x) -> !in(x, inter), p1)
    filter!((x) -> !in(x, inter), p2)
    p = [p1, reverse(p2)]
    return sum(distanceof, p)
end


"""
Find the number of edges in the shortest path between two nodes in a tree.

**Parameters:**
* `tree`: The Phylogeny to search in.
* `n1`: The first node.
* `n2`: The second node.
"""
function depth(tree::Phylogeny, n1::PhyNode, n2::PhyNode)
    p = pathbetween(tree, n1, n2)
    return length(p) == 1 ? 0 : length(p) - 1
end


"""
Find the distance between a node and the root of a tree.

**Parameters:**
* `tree`: The Phylogeny to search in.
* `n1`:   The node.
"""
function distance(tree::Phylogeny, n1::PhyNode)
    p = Tip2Root(n1)
    return sum(getbranchlength, p)
end


"""
Find the distance of each node from the root.

**Parameters:**
* `tree`: The Phylogeny to measure.
"""
function distance(tree::Phylogeny)
    distances = TreeAnnotations(tree, Float64)
    function updatedistances(node::PhyNode, currentdist::Float64)
        distances[node] = currentdist
        for child in children(node)
            updatedistances(child, currentdist + distanceof(child))
        end
    end
    updatedistances(root(tree), distanceof(root(tree)))
    return distances
end


"""
Find the depth of each node from the root.

**Parameters:**
* `tree`: The Phylogeny to measure.
"""
function depth(tree::Phylogeny)
    depths = TreeAnnotations(tree, Int)
    function updatedepths(node::PhyNode, currentdepth::Int)
        depths[node] = currentdepth
        for child in children(node)
            updatedepths(child, currentdepth + 1)
        end
    end
    updatedepths(root(tree), 0)
    return depths
end


"""
Graft a Phylogeny to the node of another tree, creates a parent-child relationship between the input node, and the root of the input phylogeny.
This function sets the root of the phylogeny object to an empty node, as the root, and so the entire structure of the tree,
has been moved to the tree containing the specified parent PhyNode.

**Parameters:**
* `parent`: The PhyNode to add the root of phylogeny too.
* `child`:  The Phylogeny for which the root is to be attached to the input parent PhyNode.
"""
function graft!(parent::PhyNode, child::Phylogeny)
    graft!(parent, root(children))
    root_unsafe!(child, PhyNode())
end


"""
Graft a Phylogeny to the node of another tree, creates a parent-child relationship between the input node, and the root of the input phylogeny.

**Parameters:**
* `parent`: The PhyNode to add the root of phylogeny too.
* `child`:  The Phylogeny for which the root is to be attached to the input parent PhyNode.
* `bl`:     Branch length connecting the parent node to the grafted phylogeny.
"""
function graft!(parent::PhyNode, child::Phylogeny, bl::Float64)
    branchlength!(root(child), bl)
    graft!(parent, child)
end


"""
Set the root field of a Phylogeny variable.

**Warning** This is different from the other `root!` methods, which rearrange the structure of a Phylogeny, rooting it based on an outgroup or midpoint.
rather, this function simply alters the root field. Generally this should not be used, except as a step in other methods. Careless use of this could result in loosing part of a tree for instance.

**Parameters:**
* `tree`: The phylogeny for which the root is to be set.
* `node`: The PhyNode that is to become the root of the tree.
"""
function root_unsafe!(tree::Phylogeny, node::PhyNode)
    tree.root = node
end


"""
Ladderize a phylogeny, i.e. sort the descendants of all nodes by their size. This creates a prettier phylogeny for plotting. The nodes of the phylogeny themselves will be ladderized, also affecting other phylogenies in memory containing the same nodes. Note that this only works for rooted trees.

**Parameters:**
* `x`: The phylogeny to ladderize
"""
function ladderize!(x::Phylogeny)
    function loc!(node::PhyNode)
        if isleaf(node)
            return 1
        end

        sizes = map(loc!, children(node))
        node.children[:] = node.children[sortperm(sizes)]
        sum(sizes) + 1
    end

    loc!(phy.root)
    nothing
end

"""
Create a ladderized copy of a phylogeny Note that this only works for rooted trees.

**Parameters:**
* `x`: The phylogeny to ladderize
"""
function ladderize(x::Phylogeny)
    ret = deepcopy(x)
    ladderize!(ret)
    ret
end
