

# An IntervalCollection is an efficiently stored and indexed set of annotated
# genomic intervals. It looks something like this.
#
#                                   ┌────┐            ┌────┐
#         Each Sequence has         │chr1│            │chr2│
#           an associated       ┌───┴────┴────┐   ┌───┴────┴────┐      ...
#            IntervalTree       │IntervalTree │   │IntervalTree │
#                               └─────────────┘   └─────────────┘
#                                      │
#                         ┌────────────┴───┬────────────┐
#                 ╔═══════▼══════╗ ╔═══════▼══════╗     ▼
#                 ║(10000, 20000)║ ║(35000, 40000)║
#                 ╚══════════════╝ ╚══════════════╝     ...
#                         │
#    Intervals            ▼
#     map to       ╔════════════╗
#   liked lists    ║ Metadata 1 ║
#   of metadata    ╚════════════╝
#                         │
#                         ▼
#                  ╔════════════╗
#                  ║ Metadata 2 ║
#                  ╚════════════╝
#                         │
#                         ▼
#
#                        ...
#
# Strand information is stored in the metadata.
#
typealias IntervalCollectionTree{T} IntervalTree{Int64, Interval{T}}

type IntervalCollection{T} <: IntervalStream{T}
    # Sequence name mapped to IntervalTree, which in turn maps intervals to
    # a list of metadata.
    trees::Dict{String, IntervalCollectionTree{T}}

    # Keep track of the number of stored intervals
    length::Int

    # A vector of values(trees) sorted on sequence name.
    # This is used to iterate intervals as efficiently as possible, but is only
    # updated as needed, indicated by the ordered_trees_outdated flag.
    ordered_trees::Vector{IntervalCollectionTree{T}}
    ordered_trees_outdated::Bool


    function IntervalCollection()
        return new(Dict{String, IntervalCollectionTree{T}}(), 0,
                   IntervalCollectionTree{T}[], false)
    end

    # TODO: bulk insertion
    #function IntervalCollection(intervals::Vector{T})
    #end
end


function update_ordered_trees!{T}(ic::IntervalCollection{T})
    if ic.ordered_trees_outdated
        ic.ordered_trees = collect(IntervalCollectionTree{T}, values(ic.trees))
        p = sortperm(collect(String, keys(ic.trees)), lt=alphanum_isless)
        ic.ordered_trees = ic.ordered_trees[p]
        ic.ordered_trees_outdated = false
    end
end


function push!{T}(ic::IntervalCollection{T}, i::Interval{T})
    if !haskey(ic.trees, i.seqname)
        tree = IntervalCollectionTree{T}()
        ic.trees[i.seqname] = tree
        ic.ordered_trees_outdated = true
    else
        tree = ic.trees[i.seqname]
    end

    push!(tree, i)
    ic.length += 1
end


function show(io::IO, ic::IntervalCollection)
    const max_entries = 8
    n_entries = length(ic)
    println(io, "IntervalCollection with $(n_entries) intervals:")
    if n_entries > 0
        for (k, i) in enumerate(ic)
            if k > 8
                break
            end
            println(io, "  ", i)
        end
        if n_entries > max_entries
            print(io, "  ⋮")
        end
    end
end


function length(ic::IntervalCollection)
    return ic.length
end


typealias IntervalCollectionTreeIteratorState{T} IntervalTrees.IntervalBTreeIteratorState{Int64, Interval{T}, 64}


immutable IntervalCollectionIteratorState{T}
    i::Int # index into ordered_trees
    tree_state::IntervalCollectionTreeIteratorState{T}

    function IntervalCollectionIteratorState(i::Int)
        return new(i)
    end

    function IntervalCollectionIteratorState(i::Int, tree_state)
        return new(i, tree_state)
    end
end


function start{T}(ic::IntervalCollection{T})
    update_ordered_trees!(ic)
    i = 1
    while i <= length(ic.ordered_trees)
        tree_state = start(ic.ordered_trees[i])
        if !done(ic.ordered_trees[i], tree_state)
            return IntervalCollectionIteratorState{T}(i, tree_state)
        end
        i += 1
    end

    return IntervalCollectionIteratorState{T}(i)
end


function next{T}(ic::IntervalCollection{T},
                      state::IntervalCollectionIteratorState{T})
    i = state.i
    value, tree_state = next(ic.ordered_trees[i], state.tree_state)

    if done(ic.ordered_trees[i], tree_state)
        i += 1
        while i <= length(ic.ordered_trees)
            tree_state = start(ic.ordered_trees[i])
            if !done(ic.ordered_trees[i], tree_state)
                break
            end
            i += 1
        end
    end

    return value, IntervalCollectionIteratorState{T}(i, tree_state)
end


function done{T}(ic::IntervalCollection{T},
                      state::IntervalCollectionIteratorState{T})
    return state.i > length(ic.ordered_trees)
end


immutable IntersectIterator{S, T}
    a_trees::Vector{IntervalCollectionTree{S}}
    b_trees::Vector{IntervalCollectionTree{T}}
end


immutable IntersectIteratorState{S, T}
    i::Int # index into a_trees/b_trees.
    intersect_iterator
    intersect_iterator_state

    function IntersectIteratorState(i::Int)
        return new(i)
    end

    function IntersectIteratorState(i::Int, intersect_iterator,
                                    intersect_iterator_state)
        return new(i, intersect_iterator, intersect_iterator_state)
    end
end


@doc """
Iterate over pairs of intersecting intervals in two IntervalCollections.
""" ->
function intersect{S, T}(a::IntervalCollection{S}, b::IntervalCollection{T})
    seqnames = collect(String, intersect(Set(keys(a.trees)), Set(keys(b.trees))))
    sort!(seqnames, lt=alphanum_isless)

    a_trees = IntervalCollectionTree{S}[a.trees[seqname] for seqname in seqnames]
    b_trees = IntervalCollectionTree{T}[b.trees[seqname] for seqname in seqnames]

    return IntersectIterator{S, T}(a_trees, b_trees)
end


function start{S, T}(it::IntersectIterator{S, T})
    i = 1
    while i <= length(it.a_trees)
        intersect_iterator = intersect(it.a_trees[i], it.b_trees[i])
        intersect_iterator_state = start(intersect_iterator)
        if !done(intersect_iterator, intersect_iterator_state)
            return IntersectIteratorState{S, T}(i, intersect_iterator,
                                                intersect_iterator_state)
        end
        i += 1
    end

    return IntersectIteratorState{S, T}(i)
end


function next{S, T}(it::IntersectIterator{S, T},
                         state::IntersectIteratorState{S, T})
    intersect_iterator = state.intersect_iterator
    value, intersect_iterator_state = next(intersect_iterator,
                                           state.intersect_iterator_state)
    i = state.i
    if done(intersect_iterator, intersect_iterator_state)
        i += 1
        while i <= length(it.a_trees)
            intersect_iterator = intersect(it.a_trees[i], it.b_trees[i])
            intersect_iterator_state = start(intersect_iterator)
            if !done(intersect_iterator, intersect_iterator_state)
                break
            end
            i += 1
        end
    end

    return value, IntersectIteratorState{S, T}(i, intersect_iterator,
                                               intersect_iterator_state)
end


function done{S, T}(it::IntersectIterator{S, T},
                         state::IntersectIteratorState{S, T})
    return state.i > length(it.a_trees)
end

