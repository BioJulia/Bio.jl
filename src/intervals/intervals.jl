
module Intervals

import Iterators
using Base.Intrinsics
using DataStructures
using Docile
using IntervalTrees

export Strand, Interval, IntervalCollection,
       STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH,
       isoverlapping

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


@doc """
A genomic interval specifies interval with some associated metadata.
""" ->
immutable Interval{T}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end


function Base.isless{T}(a::Interval{T}, b::Interval{T})
    if a.seqname != b.seqname
        return alphanum_isless(a.seqname, b.seqname)
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


# In IntervalCollection, all the data associated with an interval is stored in a linked
# list.
abstract IntervalMetadataList{T}
immutable IntervalMetadataNil{T} <: IntervalMetadataList{T} end


immutable IntervalMetadataNode{T} <: IntervalMetadataList{T}
    strand::Strand
    metadata::T
    next::IntervalMetadataList{T}
end


function Base.start(node::IntervalMetadataList)
    node
end

function Base.next{T}(::IntervalMetadataList{T}, state::IntervalMetadataList{T})
    return state, state.next
end

function Base.done{T}(::IntervalMetadataList{T}, state::IntervalMetadataList{T})
    return isa(state, IntervalMetadataNil{T})
end


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
typealias IntervalCollectionTree{T} IntervalTree{Int64, IntervalMetadataNode{T}}

type IntervalCollection{T}
    # Sequence name mapped to IntervalTree, which in turn maps intervals to
    # a list of metadata.
    trees::SortedDict{String, IntervalCollectionTree{T}, Base.Order.Lt}

    # Keep track of the number of stored intervals
    length::Int

    function IntervalCollection()
        return new(
            SortedDict(Dict{String, IntervalCollectionTree{T}}(),
                       Base.Order.Lt(alphanum_isless)), 0)
    end
end


function Base.push!{T}(is::IntervalCollection{T}, i::Interval{T})
    if !haskey(is.trees, i.seqname)
        tree = IntervalTree{Int64, IntervalMetadataNode{T}}()
        is.trees[i.seqname] = tree
    else
        tree = is.trees[i.seqname]
    end

    if !haskey(tree, (i.first, i.last))
        node = IntervalMetadataNode{T}(i.strand, i.metadata, IntervalMetadataNil{T}())
    else
        node = IntervalMetadataNode{T}(i.strand, i.metadata, tree[(i.first, i.last)])
    end

    is.length += 1
    tree[(i.first, i.last)] = node
end


function Base.show(io::IO, is::IntervalCollection)
    const max_entries = 8
    n_entries = length(is)
    println(io, "IntervalCollection with $(n_entries) intervals:")
    if n_entries > 0
        for (k, i) in enumerate(is)
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


function Base.length(is::IntervalCollection)
    return is.length
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


# Iteration over entries in IntervalCollection. Unfortunately, this is fairly
# complicated since there are several levels we have to iterate over:
# each metadata in each entry in each IntervalTree.

typealias IntervalCollectionTreeIteratorState{T} IntervalTrees.IntervalBTreeIteratorState{Int64, IntervalMetadataNode{T}, 64}

immutable IntervalCollectionIteratorState{T}
    trees_state::DataStructures.SDIterationState{String, IntervalCollectionTree{T}, Base.Order.Lt}
    interval::((Int64, Int64), IntervalMetadataList{T})
    seqname::String
    tree::IntervalCollectionTree{T}
    tree_state::IntervalCollectionTreeIteratorState{T}

    function IntervalCollectionIteratorState(
                  trees_state::DataStructures.SDIterationState{
                        String, IntervalCollectionTree{T}, Base.Order.Lt})
        return new(trees_state, ((-1, -1), IntervalMetadataNil{T}()))
    end

    function IntervalCollectionIteratorState(
                  trees_state::DataStructures.SDIterationState{
                        String, IntervalCollectionTree{T}, Base.Order.Lt},
                  interval::((Int64, Int64), IntervalMetadataList{T}),
                  seqname::String,
                  tree::IntervalCollectionTree{T},
                  tree_state::IntervalCollectionTreeIteratorState{T})
        return new(trees_state, interval, seqname, tree, tree_state)
    end
end


function Base.start{T}(ic::IntervalCollection{T})
    trees_state = start(ic.trees)
    if done(ic.trees, trees_state)
        return IntervalCollectionIteratorState{T}(trees_state)
    end

    (seqname, tree), trees_state = next(ic.trees, trees_state)
    tree_state = start(tree)
    interval, tree_state = next(tree, tree_state)
    return IntervalCollectionIteratorState{T}(trees_state, interval, seqname,
                                              tree, tree_state)
end


function Base.next{T}(ic::IntervalCollection{T}, state::IntervalCollectionIteratorState{T})
    first, last = state.interval[1]
    interval_metadata = state.interval[2]

    i = Interval{T}(state.seqname, first, last, interval_metadata.strand,
                    interval_metadata.metadata)

    trees_state = state.trees_state
    seqname     = state.seqname
    tree        = state.tree
    tree_state  = state.tree_state
    interval_metadata = interval_metadata.next
    if isa(interval_metadata, IntervalMetadataNil{T})
        if done(tree, tree_state)
            if done(ic.trees, trees_state)
                @goto end_of_iterator
            end
            (seqname, tree), trees_state = next(ic.trees, trees_state)
            tree_state = start(tree)
        end
        ((first, last), interval_metadata), tree_state = next(tree, tree_state)
    end
    @label end_of_iterator

    newstate = IntervalCollectionIteratorState{T}(trees_state,
                                                  ((first, last), interval_metadata),
                                                  seqname, tree, tree_state)

    return i, newstate
end


function Base.done{T}(ic::IntervalCollection{T}, state::IntervalCollectionIteratorState{T})
    return done(ic.trees, state.trees_state) && isa(state.interval[2], IntervalMetadataNil{T})
end


immutable IntersectIterator{S, T}
    seqnames::Vector{String}
    a::IntervalCollection{S}
    b::IntervalCollection{T}
end


@doc """
Iterate over pairs of intersecting intervals.
""" ->
function Base.intersect{S, T}(a::IntervalCollection{S}, b::IntervalCollection{T})
    seqnames = collect(String, intersect(Set(keys(a.trees)), Set(keys(b.trees))))
    sort!(seqnames, lt=alphanum_isless)

    return IntersectIterator{S, T}(seqnames, a, b)
end


immutable IntersectIteratorState{S, T}
    seqname_idx::Int
    tree_intersect_it
    tree_intersect_it_state
    metadata_product
    metadata_product_state
    first_a::Int
    last_a::Int
    first_b::Int
    last_b::Int

    function IntersectIteratorState(seqname_idx::Int)
        return new(seqname_idx)
    end

    function IntersectIteratorState(seqname_idx::Int,
                                    tree_intersect_it,
                                    tree_intersect_it_state,
                                    metadata_product,
                                    metadata_product_state,
                                    first_a, last_a,
                                    first_b, last_b)
        return new(seqname_idx,
                   tree_intersect_it, tree_intersect_it_state,
                   metadata_product, metadata_product_state,
                   first_a, last_a, first_b, last_b)
    end
end


function Base.start{S, T}(it::IntersectIterator{S, T})
    seqname_idx = 1
    if seqname_idx > length(it.seqnames)
        return IntersectIteratorState{S, T}(seqname_idx)
    end

    seqname = it.seqnames[seqname_idx]
    tree_a = it.a.trees[seqname]
    tree_b = it.b.trees[seqname]

    tree_intersect_it = intersect(tree_a, tree_b)
    tree_intersect_it_state = start(tree_intersect_it)

    while done(tree_intersect_it, tree_intersect_it_state)
        seqname_idx += 1
        if seqname_idx > length(it.seqnames)
            break
        end
        seqname = it.seqnames[seqname_idx]
        tree_a = it.a.trees[seqname]
        tree_b = it.b.trees[seqname]

        tree_intersect_it = intersect(tree_a, tree_b)
        tree_intersect_it_state = start(tree_intersect_it)
    end

    if done(tree_intersect_it, tree_intersect_it_state)
        return IntersectIteratorState{S, T}(seqname_idx)
    end

    intersect_value, tree_intersect_it_state =
        next(tree_intersect_it, tree_intersect_it_state)
    # TODO: Returning these complex tuples is probably a terrible idea
    (((first_a, last_a), metadata_list_a),
     ((first_b, last_b), metadata_list_b)) = intersect_value

    metadata_product = Iterators.product(metadata_list_a, metadata_list_b)
    metadata_product_state = start(metadata_product)

    return IntersectIteratorState{S, T}(seqname_idx,
                                        tree_intersect_it,
                                        tree_intersect_it_state,
                                        metadata_product,
                                        metadata_product_state,
                                        first_a, last_a,
                                        first_b, last_b)
end


function Base.next{S, T}(it::IntersectIterator{S, T},
                         state::IntersectIteratorState{S, T})
    (metadata_a, metadata_b), metadata_product_state =
        next(state.metadata_product, state.metadata_product_state)
    seqname = it.seqnames[state.seqname_idx]
    interval_a = Interval{S}(seqname, state.first_a, state.last_a, metadata_a.strand,
                             metadata_a.metadata)
    interval_b = Interval{T}(seqname, state.first_b, state.last_b, metadata_b.strand,
                             metadata_b.metadata)
    value = (interval_a, interval_b)

    seqname_idx = state.seqname_idx
    metadata_product = state.metadata_product
    tree_intersect_it = state.tree_intersect_it
    tree_intersect_it_state = state.tree_intersect_it_state
    first_a, last_a = state.first_a, state.last_a
    first_b, last_b = state.first_b, state.last_b

    while done(metadata_product, metadata_product_state)
        if !done(tree_intersect_it, tree_intersect_it_state)
            intersect_value, tree_intersect_it_state =
                next(tree_intersect_it, tree_intersect_it_state)

            (((first_a, last_a), metadata_list_a),
            ((first_b, last_b), metadata_list_b)) = intersect_value

            metadata_product = Iterators.product(metadata_list_a, metadata_list_b)
            metadata_product_state = start(metadata_product)
            break
        else
            seqname_idx += 1
            if seqname_idx > length(it.seqnames)
                break
            end
            seqname = it.seqnames[seqname_idx]
            tree_a = it.a.trees[seqname]
            tree_b = it.b.trees[seqname]

            tree_intersect_it = intersect(tree_a, tree_b)
            tree_intersect_it_state = start(tree_intersect_it)
        end
    end

    state = IntersectIteratorState{S, T}(
                seqname_idx,
                tree_intersect_it,
                tree_intersect_it_state,
                metadata_product,
                metadata_product_state,
                first_a, last_a,
                first_b, last_b)

    return value, state
end

function Base.done{S, T}(it::IntersectIterator{S, T},
                         state::IntersectIteratorState{S, T})
    return state.seqname_idx > length(it.seqnames)
end


end

