
module Intervals

using Base.Intrinsics
using IntervalTrees
using DataStructures

export Strand, Interval, IntervalSet, STRAND_NA, STRAND_POS, STRAND_NEG, STRAND_BOTH

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


immutable Interval{T}
    seqname::String
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end


function Base.isless{T}(a::Interval{T}, b::Interval{T})
    if a.seqname != b.seqname
        return a.seqname < b.seqname
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


function Base.show(io::IO, i::Interval)
    print(io, i.seqname, ":", i.first, "-", i.last, "    ", i.strand, "    ", i.metadata)
end


# In IntervalSet, all the data associated with an interval is stored in a linked
# list.
abstract IntervalMetadataList{T}
immutable IntervalMetadataNil{T} <: IntervalMetadataList{T} end


immutable IntervalMetadataNode{T} <: IntervalMetadataList{T}
    strand::Strand
    metadata::T
    next::IntervalMetadataList{T}
end


# An IntervelSet is an efficiently stored and indexed set of annotated genomic
# intervals. It looks something like this.
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
typealias IntervalSetTree{T} IntervalTree{Int64, IntervalMetadataNode{T}}

type IntervalSet{T}
    # Sequence name mapped to IntervalTree, which in turn maps intervals to
    # a list of metadata.
    trees::Dict{String, IntervalSetTree{T}}

    # Keep track of the number of stored intervals
    length::Int

    function IntervalSet()
        return new(Dict{String, IntervalSetTree{T}}(), 0)
    end
end


function Base.push!{T}(is::IntervalSet{T}, i::Interval{T})
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


function Base.show(io::IO, is::IntervalSet)
    const max_entries = 8
    println(io, "IntervalSet with $(length(is)) intervals:")
    for (k, i) in enumerate(is)
        if k > 8
            break
        end
        println(io, "  ", i)
    end
    print(io, "  ⋮")
end


function Base.length(is::IntervalSet)
    return is.length
end


# Iteration over entries in IntervalSet. Unfortunately, this is fairly
# complicated since there are several levels we have to iterate over:
# each metadata in each entry in each IntervalTree.

typealias IntervalSetTreeIteratorState{T} IntervalTrees.IntervalBTreeIteratorState{Int64, IntervalMetadataNode{T}, 64}

immutable IntervalSetIteratorState{T}
    trees_state::Int
    seqname::String
    tree::IntervalSetTree{T}
    tree_state::IntervalSetTreeIteratorState{T}
    interval::((Int64, Int64), IntervalMetadataList{T})

    function IntervalSetIteratorState(trees_state::Int)
        return new(trees_state)
    end

    function IntervalSetIteratorState(trees_state::Int,
                                      seqname::String,
                                      tree::IntervalSetTree{T},
                                      tree_state::IntervalSetTreeIteratorState{T},
                                      interval::((Int64, Int64), IntervalMetadataList{T}))
        return new(trees_state, seqname, tree, tree_state, interval)
    end
end


function Base.start{T}(is::IntervalSet{T})
    trees_state = start(is.trees)
    if done(is.trees, trees_state)
        return IntervalSetIteratorState{T}()
    end

    (seqname, tree), trees_state = next(is.trees, trees_state)
    tree_state = start(tree)
    interval, tree_state = next(tree, tree_state)
    return IntervalSetIteratorState{T}(trees_state, seqname, tree,
                                       tree_state, interval)
end


function Base.next{T}(is::IntervalSet{T}, state::IntervalSetIteratorState{T})
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
            if done(is.trees, trees_state)
                @goto end_of_iterator
            end
            (seqname, tree), trees_state = next(is.trees, trees_state)
            tree_state = start(tree)
        end
        ((first, last), interval_metadata), tree_state = next(tree, tree_state)
    end
    @label end_of_iterator

    newstate = IntervalSetIteratorState{T}(trees_state, seqname, tree,
                                           tree_state, ((first, last),
                                           interval_metadata))
    return i, newstate
end


function Base.done{T}(is::IntervalSet{T}, state::IntervalSetIteratorState{T})
    return isa(state.interval[2], IntervalMetadataNil{T})
end


end

