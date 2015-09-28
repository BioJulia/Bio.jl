# =================
# AlignmentAnchors
# =================


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
    alnpos::Int
    srcpos::Int
    op::Operation
    function AlignmentAnchor(gp::Int = 0, sp::Int = 0, op::Operation = OP_INVALID)
        return new(gp, sp, op)
    end
end


# Basic operators for AlignmentAnchors
# -------------------------------------

"""
Print to screen or other IO stream, a formatted description of an
AlignmentAnchor.
"""
function show(io::IO, anc::AlignmentAnchor)
    write(io, "Alignment Position: $(anc.alnpos), Source Position: $(anc.srcpos), Operation: $(anc.op)")
end


"""
Copy an AlignmentAnchor to a new AlignmentAnchor variable.

Takes only one parameter of AlignmentAnchor type.
"""
function copy(src::AlignmentAnchor)
    return AlignmentAnchor(src.alnpos, src.srcpos, src.op)
end


# Basic operators for boolean operations consider
# positions, not operations.


"""
Check for equity of two AlignmentAnchors.
"""
function ==(a::AlignmentAnchor, b::AlignmentAnchor)
    return a.alnpos == b.alnpos && a.srcpos == b.srcpos
end


"""
Check for inequity of two AlignmentAnchors.
"""
function !=(a::AlignmentAnchor, b::AlignmentAnchor)
    return !(a == b)
end


"""
Check that AlignmentAnchor a, is less than AlignmentAnchor b.
"""
function <(a::AlignmentAnchor, b::AlignmentAnchor)
    return a.alnpos < b.alnpos || a.srcpos < b.srcpos
end


"""
Check that AlignmentAnchor a, is greater than AlignmentAnchor b.
"""
function >(a::AlignmentAnchor, b::AlignmentAnchor)
    return a.alnpos > b.alnpos || a.srcpos > b.srcpos
end


"""
Check that AlignmentAnchor a, is less than or equal to AlignmentAnchor b.
"""
function <=(a::AlignmentAnchor, b::AlignmentAnchor)
    return a.alnpos <= b.alnpos || a.srcpos < b.srcpos
end


"""
Check that AlignmentAnchor a, is greater than or equal to AlignmentAnchor b.
"""
function >=(a::AlignmentAnchor, b::AlignmentAnchor)
    return a.alnpos >= b.alnpos || a.srcpos > b.srcpos
end


"""
Check whether the alignment anchor contains the specified operation.
"""
function hasop(anc::AlignmentAnchor, op::Operation)
    return anc.op == op
end



# AlignmentAnchors Definition
# ----------------------------

typealias AlignmentAnchors Vector{AlignmentAnchor}

function show(io::IO, aa::AlignmentAnchors)
    out = ""
    for i in aa
        out *= "($(i.srcpos), $(i.alnpos), $(Char(i.op)))"
    end
    write(io, out)
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

function upper_bound_anchor(arr::AlignmentAnchors, i::Int, o::AlnPosOrdering)
    return searchsortedlast(arr, i, o) + 1
end

function lower_bound_anchor(arr::AlignmentAnchors, i::Int, o::AlnPosOrdering)
    return searchsortedfirst(arr, i, o)
end


"""
Returns the indicies of the anchors that
"""
function findposition(anchors::AlignmentAnchors, position::Int, ordering::AlnPosOrdering)
    lobracket = searchsortedlast(anchors, position, ordering)
    hibracket = lobracket + 1
    # We probably need code here to handle a few edge cases, or have this as an
    # unsafe function and rely on calling functions to make sure edge cases are
    # handled.
    return lobracket, hibracket
end


immutable AlignedSequence
    src # Sequence from Bio.seq.
    anchors::AlignmentAnchors
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
