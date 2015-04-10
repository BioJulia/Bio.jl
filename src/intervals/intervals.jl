
module Intervals

import Iterators
using Base.Intrinsics
#using DataStructures
using Docile
using IntervalTrees

export Strand, Interval, IntervalCollection,
       STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH,
       isoverlapping

include("stream_buffer.jl")

bitstype 8 Strand

convert(::Type{Strand}, strand::Uint8) = box(Strand, unbox(Uint8, strand))
convert(::Type{Uint8}, strand::Strand) = box(Uint8, unbox(Strand, strand))

const STRAND_NA   = convert(Strand, 0b000)
const STRAND_POS  = convert(Strand, 0b001)
const STRAND_NEG  = convert(Strand, 0b010)
const STRAND_BOTH = convert(Strand, 0b011)

function Base.show(io::IO, strand::Strand)
    if strand == STRAND_NA
        print(io, "(indeterminate strand)")
    elseif strand == STRAND_POS
        print(io, "+")
    elseif strand == STRAND_NEG
        print(io, "-")
    elseif strand == STRAND_BOTH
        print(io, ".")
    else
        print(io, "(undefined strand)")
    end
end


function Base.isless(a::Strand, b::Strand)
    return convert(Uint8, a) < convert(Uint8, b)
end


# Note, just to be clear: this shadows IntervalTrees.Interval
@doc """
A genomic interval specifies interval with some associated metadata.
""" ->
immutable Interval{T} <: IntervalTrees.AbstractInterval{Int64}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end


function IntervalTrees.first(i::Interval)
    return i.first
end


function IntervalTrees.last(i::Interval)
    return i.last
end


function Base.isless{T}(a::Interval{T}, b::Interval{T},
                        seqname_isless::Function=alphanum_isless)
    if a.seqname != b.seqname
        return seqname_isless(a.seqname, b.seqname)
    elseif a.first != b.first
        return a.first < b.first
    elseif a.last != b.last
        return a.last < b.last
    elseif a.strand != b.strand
        return a.strand < b.strand
    else
        return a.metadata < b.metadata
    end
end


function precedes{T}(a::Interval{T}, b::Interval{T},
                     seqname_isless::Function=alphanum_isless)
    return seqname_isless(a.seqname, b.seqname) ||
           a.seqname == b.seqname && a.last < b.first
end


function =={T}(a::Interval{T}, b::Interval{T})
    return a.seqname  == b.seqname &&
           a.first    == b.first &&
           a.last     == b.last &&
           a.strand   == b.strand &&
           a.metadata == b.metadata
end


function isoverlapping{S, T}(a::Interval{S}, b::Interval{T})
    return a.seqname == b.seqname && a.first <= b.last && b.first <= a.last
end


function Base.show(io::IO, i::Interval)
    print(io, i.seqname, ":", i.first, "-", i.last, "    ", i.strand, "    ", i.metadata)
end


@doc """
A type deriving `IntervalStream{T}` must be iterable and produce
Interval{T} objects in sorted order.
""" ->
abstract IntervalStream{T}


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


function Base.length(ic::IntervalCollection)
    return ic.length
end


@doc """
A comparison function used to sort on numbers within text.

This is useful since sequences are often named things like "chr12" or
"read1234". Treating the numbers as numbers and not text gives a more natural
ordering.

This is similar to the '--version-sort' option in GNU coreutils sort.
""" ->
function alphanum_isless(a::String, b::String)
    i = 1
    j = 1

    # match up to the first digit
    k0 = 0 # position of first digit
    while i <= length(a) && j <= length(b)
        if isdigit(a[i]) && isdigit(b[j])
            k0 = i
            break
        else
            if a[i] != b[j]
                return a[i] < b[j]
            end
        end
        i = nextind(a, i)
        j = nextind(b, j)
    end

    # match numbers
    ka1, kb1 = 0, 0
    while i <= length(a) && isdigit(a[i])
        ka1 = i
        i = nextind(a, i)
    end
    while j <= length(b) && isdigit(b[j])
        kb1 = j
        j = nextind(b, j)
    end

    if ka1 == 0 && kb1 != 0
        return true
    elseif ka1 != 0 && kb1 == 0
        return false
    elseif ka1 != 0 && kb1 != 0
        aval = parse(Int, a[k0:ka1])
        bval = parse(Int, b[k0:kb1])
        if aval != bval
            return aval < bval
        end
    end

    # match suffixes
    while i <= length(a) && j <= length(b)
        if a[i] != b[j]
            return a[i] < b[j]
        end
        i = nextind(a, i)
        j = nextind(b, j)
    end

    return j <= length(b)
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
function Base.intersect{S, T}(a::IntervalCollection{S}, b::IntervalCollection{T})
    seqnames = collect(String, intersect(Set(keys(a.trees)), Set(keys(b.trees))))
    sort!(seqnames, lt=alphanum_isless)

    a_trees = IntervalCollectionTree{S}[a.trees[seqname] for seqname in seqnames]
    b_trees = IntervalCollectionTree{T}[b.trees[seqname] for seqname in seqnames]

    return IntersectIterator{S, T}(a_trees, b_trees)
end


function Base.start{S, T}(it::IntersectIterator{S, T})
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


function Base.next{S, T}(it::IntersectIterator{S, T},
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


function Base.done{S, T}(it::IntersectIterator{S, T},
                         state::IntersectIteratorState{S, T})
    return state.i > length(it.a_trees)
end


# IntervalStreams
# ---------------
#
# Often we'd rather avoid reading a large interval dataset into an
# IntervalCollection. We can still do efficient intersections if we can assume
# the data is sorted


type IntervalStreamIntersectIterator{S, T}
    a::IntervalStream{S}
    b::IntervalStream{T}
    seqname_isless::Function

    # buffered intervals from a
    a_buffer::StreamBuffer{Interval{S}}
    b_interval::Interval{T}

    function IntervalStreamIntersectIterator(a::IntervalStream{S},
                                             b::IntervalStream{T},
                                             seqname_isless::Function)
        return new(a, b, seqname_isless, StreamBuffer{Interval{S}}())
    end
end


immutable IntervalStreamIntersectIteratorState{U, V}
    a_buffer_pos::Int
    a_state::U
    b_state::V
end


@doc """
Intersect two `IntervalStreams` returning an iterator over all pairs of
intersecting intervals.
""" ->
function Base.intersect{S, T}(a::IntervalStream{S}, b::IntervalStream{T},
                              seqname_isless::Function=isless)
    return IntervalStreamIntersectIterator{S, T}(a, b, seqname_isless)
end


function first_intersection{S, T, U, V}(a::IntervalStream{S},
                                        b::IntervalStream{T},
                                        a_state::U,
                                        b_state::V,
                                        a_interval::Interval{S},
                                        b_interval::Interval{T},
                                        seqname_isless::Function)

    while !done(a, a_state) || !done(b, b_state)
        if precedes(a_interval, b_interval) && !done(a, a_state)
            a_interval_next, a_state = next(a, a_state)
            if isless(a_interval_next, a_interval, seqname_isless)
                error("Interval streams must be sorted to perform intersection.")
            end
            a_interval = a_interval_next
        elseif precedes(b_interval, a_interval) && !done(b, b_state)
            b_interval_next, b_state = next(b, b_state)
            if isless(b_interval_next, b_interval, seqname_isless)
                error("Interval streams must be sorted to perform intersection.")
            end
            b_interval = b_interval_next
        else # intersection
            break
        end
    end

    return a_interval, b_interval, a_state, b_state
end


function Base.start{S, T}(it::IntervalStreamIntersectIterator{S, T})
    a_state = start(it.a)
    b_state = start(it.b)

    if done(it.a, a_state) || done(it.b, b_state)
        return IntervalStreamIntersectIteratorState(a_state, b_state)
    end

    a_interval, a_state = next(it.a, a_state)
    b_interval, b_state = next(it.b, b_state)
    a_interval, b_interval, a_state, b_state = first_intersection(
        it.a, it.b, a_state, b_state, a_interval, b_interval, it.seqname_isless)

    if isoverlapping(a_interval, b_interval)
        push!(it.a_buffer, a_interval)
        it.b_interval = b_interval
    end

    return IntervalStreamIntersectIteratorState(1, a_state, b_state)
end


# Do we need a state if this is going to be stateful. State can maybe just
# be the current interval in b
function Base.next{S, T, U, V}(it::IntervalStreamIntersectIterator{S, T},
                               state::IntervalStreamIntersectIteratorState{U, V})
    value = (it.a_buffer[state.a_buffer_pos], it.b_interval)
    a_state = state.a_state
    b_state = state.b_state

    # another intersection in a_buffer?
    a_buffer_pos = state.a_buffer_pos + 1
    while a_buffer_pos <= length(it.a_buffer) &&
          !isoverlapping(it.a_buffer[a_buffer_pos], it.b_interval)
        a_buffer_pos += 1
    end

    # look for an intersection for b_interval by reading more entries from a
    b_interval = it.b_interval
    if a_buffer_pos > length(it.a_buffer) && !done(it.a, a_state)
        while true
            a_interval, a_state = next(it.a, a_state)
            if isless(a_interval, it.a_buffer[end], it.seqname_isless)
                error("Interval streams must be sorted to perform intersection.")
            end
            push!(it.a_buffer, a_interval)

            if isoverlapping(b_interval, it.a_buffer[end]) ||
               precedes(b_interval, it.a_buffer[end]) ||
               done(it.a, a_state)
                break
            end
            a_buffer_pos += 1
        end

        if !isoverlapping(it.a_buffer[a_buffer_pos], b_interval)
            a_buffer_pos += 1
        end
    end

    # no more intersections against b_interval, try next interval from b
    while (a_buffer_pos > length(it.a_buffer) ||
           !isoverlapping(it.a_buffer[a_buffer_pos], b_interval)) && !done(it.b, b_state)

        b_interval_next, b_state = next(it.b, b_state)
        if isless(b_interval_next, b_interval, it.seqname_isless)
            error("Interval streams must be sorted to perform intersection.")
        end
        b_interval = b_interval_next

        # evict as much as possible from a_buffer
        a_interval = it.a_buffer[1]
        while !isempty(it.a_buffer) && precedes(it.a_buffer[1], b_interval, it.seqname_isless)
            a_interval = shift!(it.a_buffer)
        end

        if isempty(it.a_buffer)
            a_interval, b_interval, a_state, b_state = first_intersection(
                it.a, it.b, a_state, b_state, a_interval, b_interval, it.seqname_isless)
            it.b_interval = b_interval
            if isoverlapping(a_interval, b_interval)
                push!(it.a_buffer, a_interval)
                a_buffer_pos = 1
            end
            break
        elseif isoverlapping(it.a_buffer[1], b_interval)
            a_buffer_pos = 1
            it.b_interval = b_interval
        end
    end

    return value, IntervalStreamIntersectIteratorState{U, V}(
                        a_buffer_pos, a_state, b_state)
end


function Base.done{S, T, U, V}(it::IntervalStreamIntersectIterator{S, T},
                               state::IntervalStreamIntersectIteratorState{U, V})
    return state.a_buffer_pos > length(it.a_buffer)
end


# TODO: IntervalCollection(stream::IntervalStream) constructor.

end

