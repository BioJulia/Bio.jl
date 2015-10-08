
# Anchor definition
# ------------------
"""
A type to store the operation enocded in an alignment (CIGAR Operations).
Also stores the position in the alignment view of the sequences, and the
corresponding position in the unaltered source sequence (nucleotide or protein).

It stores the position values as Int types and the alignment operation is stored
as a type `Operation`, these are defined in the file `operations.jl` and loaded
into Bio.Align namespace as a series of global constants.

Calling the AlignmentAnchor default constructor with default arguments sets both
integer fields to 0, and the operation field to the global constant OP_INVALID.
"""
immutable AlignmentAnchor
    seqpos::Int
    refpos::Int
    op::Operation
end


# Basic operators for AlignmentAnchors
# -------------------------------------

"""
Print to screen or other IO stream, a formatted description of an
AlignmentAnchor.
"""
function show(io::IO, anc::AlignmentAnchor)
    print(io, "AlignmentAnchor(", anc.seqpos, ", ", anc.refpos, ", '", anc.op, "')")
end


# Alignment
# ---------

immutable Alignment
    anchors::Vector{AlignmentAnchor}
    firstref::Int
    lastref::Int

    function Alignment(anchors::Vector{AlignmentAnchor}, check::Bool=true)
        # optionally check coherence of the anchors
        if check
            # empty alignments are valid, representing an unaligned sequence
            if !isempty(anchors)
                if anchors[1].op != OP_START
                    error("Alignments must begin with on OP_START anchor.")
                end

                for i in 2:length(anchors)
                    if anchors[i].refpos < anchors[i-1].refpos ||
                       anchors[i].seqpos < anchors[i-1].seqpos
                        error("Alignment anchors must be sorted.")
                    end

                    op = anchors[i].op
                    # reference skip/delete operations
                    if isdeleteop(op)
                        if anchors[i].seqpos != anchors[i-1].seqpos
                            error("Invalid anchor positions for reference deletion.")
                        end
                    # reference insertion operations
                    elseif isinsertop(op)
                        if anchors[i].refpos != anchors[i-1].refpos
                            error("Invalid anchor positions for reference insertion.")
                        end
                    # match operations
                    elseif ismatchop(op)
                        if anchors[i].refpos - anchors[i-1].refpos !=
                           anchors[i].seqpos - anchors[i-1].seqpos
                            error("Invalid anchor positions for match operation.")
                        end
                    end
                end
            end
        end

        # compute first and last aligned reference positions
        firstref = 0
        for i in 1:length(anchors)
            if ismatchop(anchors[i].op)
                firstref = anchors[i-1].refpos
            end
        end

        lastref = 0
        for i in length(anchors):-1:1
            if ismatchop(anchors[i].op)
                lastref = anchors[i].refpos
            end
        end

        return new(anchors, firstref, lastref)
    end
end


function Alignment(refpos::Int, cigar::AbstractString)
    # TODO: conversion from start position and cigar string
end


function show(io::IO, aln::Alignment)
    # print a representation of the reference sequence
    anchors = aln.anchors
    for i in 2:length(anchors)
        if ismatchop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '.')
            end
        elseif isinsertop(anchors[i].op)
            for _ in anchors[i-1].seqpos+1:anchors[i].seqpos
                write(io, '-')
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '.')
            end
        end
    end
    write(io, '\n')

    # print a representation of the aligned sequence
    for i in 2:length(anchors)
        if ismatchop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                write(io, '.')
            end
        elseif isinsertop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                write(io, '.')
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '-')
            end
        end
    end
end



# Sorting AlignmentAnchors
# ------------------------

#=
In order to make use of Julia's sorting API and avoid having to write our own,
whilst ALSO avoiding slowing down things by passing functions about as arguments
to other functions, we make use of how the sort API constructs its searching and
sorting algorithms. It uses types that inherit from the Ordering abstract type,
and this type, combined with use of different Base.Order.lt methods - define how
the values in the array are compared to the query value and in what order.
So we take advantage of this and define our own types inheriting from
Base.Order.Ordering and create our own lt methods. These will then be used by
the core Julia sorting API efficiently.
=#

# Custom types inheriting from Base.Ordering

"""
An immutable field-less type, inheriting from Ordering. This is used with
Julia's sorting API functions to determine how an array of AlignmentAnchors is
sorted.

Specifically, it specifies that the sorting API functions use less-than
methods functions that compare the srcpos fields of the AlignmentAnchors.
"""
immutable SrcPosOrdering <: Ordering end

"""
An immutable field-less type, inheriting from Ordering. This is used with
Julia's sorting API functions to determine how an array of AlignmentAnchors is
sorted.

Specifically, it specifies that the sorting API functions use less-than
methods functions that compare the alnpos fields of the AlignmentAnchors.
"""
immutable AlnPosOrdering <: Ordering end

const BY_SRC = SrcPosOrdering()
const BY_ALN = AlnPosOrdering()

# Is AlignmentAnchor a, less than AlignmentAnchor b, according to the source
# sequence co-ordinate.


"""
Returns true if the srcpos field of AlignmentAnchor a, is less than the srcpos
field of AlignmentAnchor b.
"""
function lt(o::SrcPosOrdering, a::AlignmentAnchor, b::AlignmentAnchor)
    return a.srcpos < b.srcpos
end

"""
Returns true if the alnpos field of AlignmentAnchor a, is less than the alnpos
field of AlignmentAnchor b.
"""
function lt(o::AlnPosOrdering, a::AlignmentAnchor, b::AlignmentAnchor)
    return a.alnpos < b.alnpos
end

"""
Returns true if the srcpos field of AlignmentAnchor a, is less than the Integer
b.
"""
function lt(o::SrcPosOrdering, a::AlignmentAnchor, b::Int)
    return a.srcpos < b
end


"""
Returns true if the alnpos field of AlignmentAnchor a, is less than the Integer
b.
"""
function lt(o::AlnPosOrdering, a::AlignmentAnchor, b::Int)
    return a.alnpos < b
end


"""
Returns true if the srcpos field of Integer a, is less than the srcpos field of
the AlignmentAnchor b.
"""
function lt(o::SrcPosOrdering, a::Int, b::AlignmentAnchor)
    return a < b.srcpos
end

"""
Returns true if the srcpos field of Integer a, is less than the alnpos field of
the AlignmentAnchor b.
"""
function lt(o::AlnPosOrdering, a::Int, b::AlignmentAnchor)
    return a < b.alnpos
end

function upper_bound_anchor(arr::Vector{AlignmentAnchor}, i::Int, o::AlnPosOrdering)
    return searchsortedlast(arr, i, o) + 1
end

function lower_bound_anchor(arr::Vector{AlignmentAnchor}, i::Int, o::AlnPosOrdering)
    return searchsortedfirst(arr, i, o)
end


"""
Returns the indicies of the anchors that
"""
function findposition(anchors::Vector{AlignmentAnchor}, position::Int, ordering::AlnPosOrdering)
    lobracket = searchsortedlast(anchors, position, ordering)
    hibracket = lobracket + 1
    # We probably need code here to handle a few edge cases, or have this as an
    # unsafe function and rely on calling functions to make sure edge cases are
    # handled.
    return lobracket, hibracket
end


immutable AlignedSequence{S <: Sequence}
    seq::S
    aln::Alignment
end


function AlignedSequence{S <: Sequence}(seq::S, aln::Alignment)
    return AlignedSequence{S}(seq, aln)
end


function AlignedSequence{S <: Sequence}(seq::S, anchors::Vector{AlignmentAnchor},
                                        check::Bool=true)
    return AlignedSequence(seq, Alignment(anchors, check))
end


"""
First position in the reference sequence.
"""
function first(alnseq::AlignedSequence)
    return alnseq.aln.firstref
end


"""
Last position in the reference sequence.
"""
function last(alnseq::AlignedSequence)
    return alnseq.aln.lastref
end


# simple letters and dashes representation of an alignment
function show(io::IO, alnseq::AlignedSequence)
    # print a representation of the reference sequence
    anchors = alnseq.aln.anchors
    for i in 2:length(anchors)
        if ismatchop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '.')
            end
        elseif isinsertop(anchors[i].op)
            for _ in anchors[i-1].seqpos+1:anchors[i].seqpos
                write(io, '-')
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '.')
            end
        end
    end
    write(io, '\n')

    for i in 2:length(anchors)
        if ismatchop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                print(io, alnseq.seq[i])
            end
        elseif isinsertop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                print(io, alnseq.seq[i])
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '-')
            end
        end
    end
end


function cigar(anchors::Vector{AlignmentAnchor})
    # TODO: encode an alignment as a cigar string
end


function alntosrc(aligned_seq::AlignedSequence, alnposition::Int)
    lowidx, highidx = findposition(aligned_seq, alnposition, BY_ALN)
    high_anchor = aligned_seq.anchors[highidx]
    low_anchor = aligned_seq.anchors[lowidx]
    source_distance = high_anchor.srcpos - low_anchor.srcpos
    align_distance = alnposition - low_anchor.alnpos
    if source_distance > align_distance
        source_position = low_anchor.srcpos + align_distance
    else
        source_position = high_anchor.srcpos
    end
    return source_position
end


function srctoaln(aligned_seq::AlignedSequence, alnposition::Int)
    lowidx, highidx = findposition(aligned_seq, alnposition, BY_SRC)
end
