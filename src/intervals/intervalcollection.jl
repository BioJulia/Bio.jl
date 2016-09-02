# An IntervalCollection is an efficiently stored and indexed set of annotated
# genomic intervals. It looks something like this.
#
#                                      ┌─────┐
#                                      │trees│
#                                      └─────┘
#                                         │
#                              ┌──────────┴──────────┬────────────┐
#                              ▼                     ▼            │
#  Each sequence has       ┌──────┐              ┌──────┐         ▼
#    an associated         │ chr1 │              │ chr2 │
#    IntervalTree     ┌────┴──────┴────┐    ┌────┴──────┴────┐    ...
#                     │ IntervalTree 1 │    │ IntervalTree 2 │
#                     └────────────────┘    └────────────────┘
#                              │
#                    ┌─────────┴─────────────────────────┬────────────────────────┐
#                    │                                   │                        │
#                    ▼                                   ▼                        ▼
#   ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓  ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
#   ┃ Interval{T}(10000, 20000, ...) ┃  ┃ Interval{T}(35000, 40000, ...) ┃      ...
#   ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛  ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
#
#
#
#
#                               ┌────────────────┐
#       ordered_trees holds     │ IntervalTree 1 │
#         an array of the       ├────────────────┤
#       some IntervalTrees,     │ IntervalTree 2 │
#           ordered by          └────────────────┘
#       chromosome for fast
#       ordered iteration.             ...
#

typealias IntervalCollectionTree{T} IntervalTree{Int64,Interval{T}}

type IntervalCollection{T} <: IntervalStream{T}
    # Sequence name mapped to IntervalTree, which in turn maps intervals to
    # a list of metadata.
    trees::Dict{StringField, IntervalCollectionTree{T}}

    # Keep track of the number of stored intervals
    length::Int

    # A vector of values(trees) sorted on sequence name.
    # This is used to iterate intervals as efficiently as possible, but is only
    # updated as needed, indicated by the ordered_trees_outdated flag.
    ordered_trees::Vector{IntervalCollectionTree{T}}
    ordered_trees_outdated::Bool

    function IntervalCollection()
        return new(Dict{StringField, IntervalCollectionTree{T}}(), 0,
                   IntervalCollectionTree{T}[], false)
    end

    # bulk insertion
    function IntervalCollection(intervals::AbstractVector{Interval{T}}, sort=false)
        if sort
            sort!(intervals)
        else
            if !issorted(intervals)
                error("Intervals must be sorted, or `sort=true` set, to construct an IntervalCollection")
            end
        end

        n = length(intervals)
        trees = Dict{StringField, IntervalCollectionTree{T}}()
        i = 1
        while i <= n
            j = i
            while j <= n && intervals[i].seqname == intervals[j].seqname
                j += 1
            end
            trees[intervals[i].seqname] = IntervalCollectionTree{T}(view(intervals, i:j-1))
            i = j
        end
        return new(trees, n, IntervalCollectionTree{T}[], true)
    end
end

function IntervalCollection{T}(intervals::AbstractVector{Interval{T}}, sort=false)
    return IntervalCollection{T}(intervals, sort)
end

function IntervalCollection{T}(interval_stream::IntervalStream{T})
    intervals = collect(Interval{T}, interval_stream)
    return IntervalCollection{T}(intervals, true)
end

function IntervalCollection(intervals)
    return IntervalCollection(collect(intervals), true)
end

function update_ordered_trees!{T}(ic::IntervalCollection{T})
    if ic.ordered_trees_outdated
        ic.ordered_trees = collect(IntervalCollectionTree{T}, values(ic.trees))
        p = sortperm(collect(AbstractString, keys(ic.trees)), lt=alphanum_isless)
        ic.ordered_trees = ic.ordered_trees[p]
        ic.ordered_trees_outdated = false
    end
end

function Base.push!{T}(ic::IntervalCollection{T}, i::Interval{T})
    if !haskey(ic.trees, i.seqname)
        tree = IntervalCollectionTree{T}()
        ic.trees[i.seqname] = tree
        ic.ordered_trees_outdated = true
    else
        tree = ic.trees[i.seqname]
    end
    push!(tree, i)
    ic.length += 1
    return ic
end

function Base.show(io::IO, ic::IntervalCollection)
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

Base.length(ic::IntervalCollection) = ic.length


# Iterators
# --------

typealias IntervalCollectionTreeIteratorState{T} IntervalTrees.IntervalBTreeIteratorState{Int64,Interval{T},64}

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

function Base.start{T}(ic::IntervalCollection{T})
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

function Base.next{T}(ic::IntervalCollection{T},
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

function Base.done{T}(ic::IntervalCollection{T},
                      state::IntervalCollectionIteratorState{T})
    return state.i > length(ic.ordered_trees)
end

"""
Return an iterator over all intervals in the collection that overlap the query interval.
"""
function Base.intersect{T}(a::IntervalCollection{T}, b::Interval)
    if haskey(a.trees, b.seqname)
        return intersect(a.trees[b.seqname], b)
    else
        return IntervalTrees.IntervalIntersectionIterator{Int64, Interval{T}, 64}()
    end
end

type IntersectIterator{S, T}
    a_trees::Vector{IntervalTrees.IntervalBTree{Int64, Interval{S}, 64}}
    b_trees::Vector{IntervalTrees.IntervalBTree{Int64, Interval{T}, 64}}

    i::Int # index into a_trees/b_trees.
    intersect_iterator::IntervalTrees.IntersectionIterator{Int64, Interval{S}, 64, Interval{T}, 64}

    function IntersectIterator(a_trees, b_trees)
        return new(a_trees, b_trees)
    end
end

Base.iteratorsize(::IntersectIterator) = Base.SizeUnknown()

"Iterate over pairs of intersecting intervals in two IntervalCollections"
function Base.intersect{S, T}(a::IntervalCollection{S}, b::IntervalCollection{T})
    seqnames = collect(AbstractString, keys(a.trees) ∩ keys(b.trees))
    sort!(seqnames, lt=alphanum_isless)

    a_trees = IntervalCollectionTree{S}[a.trees[seqname] for seqname in seqnames]
    b_trees = IntervalCollectionTree{T}[b.trees[seqname] for seqname in seqnames]

    return IntersectIterator{S, T}(a_trees, b_trees)
end

function Base.start{S,T}(it::IntersectIterator{S,T})
    i = 1
    while i <= length(it.a_trees)
        intersect_iterator = intersect(it.a_trees[i], it.b_trees[i])
        intersect_iterator_state = start(intersect_iterator)
        if !done(intersect_iterator, intersect_iterator_state)
            it.i = i
            it.intersect_iterator = intersect_iterator
            return nothing
        end
        i += 1
    end
    it.i = i

    return nothing
end

function Base.next{S,T}(it::IntersectIterator{S, T}, ::Void)
    intersect_iterator = it.intersect_iterator
    value, intersect_iterator_state = next(intersect_iterator, nothing)
    i = it.i
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
    it.i = i
    it.intersect_iterator = intersect_iterator

    return value, nothing
end

function Base.done{S, T}(it::IntersectIterator{S, T}, state)
    return it.i > length(it.a_trees)
end

"""
Iterate over pairs of intersections in an IntervalCollection versus an
IntervalStream.
"""
function Base.intersect(a::IntervalCollection, b::IntervalStreamOrArray)
    return IntervalCollectionStreamIterator(metadatatype(b), a, b)
end

immutable IntervalCollectionStreamIterator{S, T, TV}
    a::IntervalCollection{S}
    b::T
end

function IntervalCollectionStreamIterator{S, T}(TV::Type, a::IntervalCollection{S},
                                                b::T)
    return IntervalCollectionStreamIterator{S, T, TV}(a, b)
end

immutable IntervalCollectionStreamIteratorState{S, TS, TV}
    intersection::IntervalTrees.Intersection{Int64, Interval{S}, 64}
    b_state::TS
    b_value::Interval{TV}

    function IntervalCollectionStreamIteratorState(intersection, b_state, b_value)
        return new(intersection, b_state, b_value)
    end

    function IntervalCollectionStreamIteratorState()
        return new(IntervalTrees.Intersection{Int64, Interval{S}, 64}())
    end
end

function Base.:(==){T}(a::IntervalCollection{T}, b::IntervalCollection{T})
    if length(a) != length(b)
        return false
    end

    for (i, j) in zip(a, b)
        if i != j
            return false
        end
    end

    return true
end

# This mostly follows from SuccessiveTreeIntersectionIterator in IntervalTrees
function Base.start{S,T,TV}(it::IntervalCollectionStreamIterator{S,T,TV})
    b_state = start(it.b)

    # TODO: We need to figure out a way to know these types at compile time.
    TS = typeof(b_state)

    intersection = IntervalTrees.Intersection{Int64, Interval{S}, 64}()
    while !done(it.b, b_state)
        b_value, b_state = next(it.b, b_state)
        if haskey(it.a.trees, b_value.seqname)
            tree = it.a.trees[b_value.seqname]
            IntervalTrees.firstintersection!(tree, b_value, intersection)
            if intersection.index != 0
                return IntervalCollectionStreamIteratorState{S, TS, TV}(
                    intersection, b_state, b_value)
            end
        end
    end

    return IntervalCollectionStreamIteratorState{S, TS, TV}()
end

function Base.next{S,T,TS,TV}(it::IntervalCollectionStreamIterator{S,T,TV},
                              state::IntervalCollectionStreamIteratorState{S,TS,TV})
    intersection = state.intersection
    entry = intersection.node.entries[intersection.index]
    return_value = (entry, state.b_value)
    b_state = state.b_state
    b_value = state.b_value
    IntervalTrees.nextintersection!(
        intersection.node, intersection.index, b_value, intersection)
    while intersection.index == 0 && !done(it.b, b_state)
        b_value, b_state = next(it.b, b_state)
        if haskey(it.a.trees, b_value.seqname)
            tree = it.a.trees[b_value.seqname]
            IntervalTrees.firstintersection!(tree, b_value, intersection)
        end
    end

    return return_value, IntervalCollectionStreamIteratorState{S, TS, TV}(
                            intersection, b_state, b_value)
end

function Base.done{S,T,TS,TV}(it::IntervalCollectionStreamIterator{S,T,TV},
                              state::IntervalCollectionStreamIteratorState{S,TS,TV})
    return state.intersection.index == 0
end
