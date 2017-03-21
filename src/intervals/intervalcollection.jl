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

# Aliases for types of IntervalTrees.jl (IC: Interval Collection).
typealias ICTree{T}                             IntervalTrees.IntervalBTree{Int64,Interval{T},64}
typealias ICTreeIteratorState{T}                IntervalTrees.IntervalBTreeIteratorState{Int64,Interval{T},64}
typealias ICTreeIntersection{T}                 IntervalTrees.Intersection{Int64,Interval{T},64}
typealias ICTreeIntersectionIterator{S,T}       IntervalTrees.IntersectionIterator{Int64,Interval{S},64,Interval{T},64}
typealias ICTreeIntervalIntersectionIterator{T} IntervalTrees.IntervalIntersectionIterator{Int64,Interval{T},64}

type IntervalCollection{T}
    # Sequence name mapped to IntervalTree, which in turn maps intervals to
    # a list of metadata.
    trees::Dict{StringField, ICTree{T}}

    # Keep track of the number of stored intervals
    length::Int

    # A vector of values(trees) sorted on sequence name.
    # This is used to iterate intervals as efficiently as possible, but is only
    # updated as needed, indicated by the ordered_trees_outdated flag.
    ordered_trees::Vector{ICTree{T}}
    ordered_trees_outdated::Bool

    function IntervalCollection()
        return new(Dict{StringField,ICTree{T}}(), 0, ICTree{T}[], false)
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
        trees = Dict{StringField, ICTree{T}}()
        i = 1
        while i <= n
            j = i
            while j <= n && intervals[i].seqname == intervals[j].seqname
                j += 1
            end
            trees[intervals[i].seqname] = ICTree{T}(view(intervals, i:j-1))
            i = j
        end
        return new(trees, n, ICTree{T}[], true)
    end
end

function IntervalCollection{T}(intervals::AbstractVector{Interval{T}}, sort=false)
    return IntervalCollection{T}(intervals, sort)
end

function IntervalCollection(intervals)
    return IntervalCollection(collect(Interval{metadatatype(intervals)}, intervals), true)
end

function update_ordered_trees!{T}(ic::IntervalCollection{T})
    if ic.ordered_trees_outdated
        ic.ordered_trees = collect(ICTree{T}, values(ic.trees))
        p = sortperm(collect(AbstractString, keys(ic.trees)), lt=isless)
        ic.ordered_trees = ic.ordered_trees[p]
        ic.ordered_trees_outdated = false
    end
end

function Base.push!{T}(ic::IntervalCollection{T}, i::Interval{T})
    if !haskey(ic.trees, i.seqname)
        tree = ICTree{T}()
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
    n_entries = length(ic)
    println(io, "IntervalCollection with $(n_entries) intervals:")
    if n_entries > 0
        for (k, i) in enumerate(ic)
            if k > 8
                break
            end
            println(io, "  ", i)
        end
        if n_entries > 8
            print(io, "  ⋮")
        end
    end
end

function Base.length(ic::IntervalCollection)
    return ic.length
end

function Base.eltype{T}(::Type{IntervalCollection{T}})
    return Interval{T}
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


# Iterators
# ---------

type IntervalCollectionIteratorState{T}
    i::Int # index into ordered_trees
    tree_state::ICTreeIteratorState{T}

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

function Base.next(ic::IntervalCollection, state)
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
    state.i, state.tree_state = i, tree_state
    return value, state
end

function Base.done(ic::IntervalCollection, state)
    return state.i > length(ic.ordered_trees)
end


# Overlaps
# --------

function eachoverlap{T}(a::IntervalCollection{T}, b::Interval)
    if haskey(a.trees, b.seqname)
        return intersect(a.trees[b.seqname], b)
    else
        return ICTreeIntervalIntersectionIterator{T}()
    end
end

function eachoverlap(a::IntervalCollection, b::IntervalCollection)
    seqnames = collect(AbstractString, keys(a.trees) ∩ keys(b.trees))
    sort!(seqnames, lt=isless)
    a_trees = [a.trees[seqname] for seqname in seqnames]
    b_trees = [b.trees[seqname] for seqname in seqnames]
    return IntersectIterator(a_trees, b_trees)
end

immutable IntersectIterator{S, T}
    a_trees::Vector{ICTree{S}}
    b_trees::Vector{ICTree{T}}
end

type IntersectIteratorState{S,T}
    i::Int  # index into a_trees/b_trees.
    intersect_iterator::ICTreeIntersectionIterator{S,T}

    function IntersectIteratorState(i)
        return new(i)
    end

    function IntersectIteratorState(i, iter)
        return new(i, iter)
    end
end

function Base.eltype{S,T}(::Type{IntersectIterator{S,T}})
    return Tuple{Interval{S},Interval{T}}
end

function Base.iteratorsize{S,T}(::Type{IntersectIterator{S,T}})
    return Base.SizeUnknown()
end

function Base.start{S,T}(it::IntersectIterator{S,T})
    i = 1
    while i <= length(it.a_trees)
        intersect_iterator = intersect(it.a_trees[i], it.b_trees[i])
        intersect_iterator_state = start(intersect_iterator)
        if !done(intersect_iterator, intersect_iterator_state)
            return IntersectIteratorState{S,T}(i, intersect_iterator)
        end
        i += 1
    end
    return IntersectIteratorState{S,T}(i)
end

function Base.next{S,T}(it::IntersectIterator{S, T}, state)
    i, intersect_iterator = state.i, state.intersect_iterator
    value, intersect_iterator_state = next(intersect_iterator, nothing)
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
    state.i, state.intersect_iterator = i, intersect_iterator
    return value, state
end

function Base.done{S, T}(it::IntersectIterator{S, T}, state)
    return state.i > length(it.a_trees)
end

function eachoverlap(a, b::IntervalCollection)
    return IntervalCollectionStreamIterator(a, b)
end

immutable IntervalCollectionStreamIterator{S,T}
    a::S
    b::IntervalCollection{T}
end

function Base.eltype{S,T}(::Type{IntervalCollectionStreamIterator{S,T}})
    return Tuple{Interval{metadatatype(S)},Interval{T}}
end

function Base.iteratorsize{S,T}(::Type{IntervalCollectionStreamIterator{S,T}})
    return Base.SizeUnknown()
end

type IntervalCollectionStreamIteratorState{Ta,Tb,U}
    intersection::ICTreeIntersection{Tb}
    a_value::Interval{Ta}
    a_state::U

    function IntervalCollectionStreamIteratorState(intersection, a_value, a_state)
        return new(intersection, a_value, a_state)
    end

    function IntervalCollectionStreamIteratorState()
        return new(ICTreeIntersection{Tb}())
    end
end

# This mostly follows from SuccessiveTreeIntersectionIterator in IntervalTrees
function Base.start{S,T}(it::IntervalCollectionStreamIterator{S,T})
    a_state = start(it.a)
    intersection = ICTreeIntersection{T}()
    while !done(it.a, a_state)
        a_value, a_state = next(it.a, a_state)
        if haskey(it.b.trees, a_value.seqname)
            tree = it.b.trees[a_value.seqname]
            IntervalTrees.firstintersection!(tree, a_value, Nullable{Interval{T}}(), intersection)
            if intersection.index != 0
                return IntervalCollectionStreamIteratorState{T,metadatatype(it.a),typeof(a_state)}(intersection, a_value, a_state)
            end
        end
    end
    return IntervalCollectionStreamIteratorState{S,metadatatype(it.a),typeof(a_state)}()
end

function Base.next{S,T}(it::IntervalCollectionStreamIterator{S,T}, state)
    intersection = state.intersection
    entry = intersection.node.entries[intersection.index]
    return_value = (state.a_value, entry)
    IntervalTrees.nextintersection!(intersection.node, intersection.index, state.a_value, intersection)
    while intersection.index == 0 && !done(it.a, state.a_state)
        state.a_value, state.a_state = next(it.a, state.a_state)
        if haskey(it.b.trees, state.a_value.seqname)
            tree = it.b.trees[state.a_value.seqname]
            IntervalTrees.firstintersection!(tree, state.a_value, Nullable{Interval{T}}(), intersection)
        end
    end
    return return_value, state
end

function Base.done(it::IntervalCollectionStreamIterator, state)
    return state.intersection.index == 0
end
