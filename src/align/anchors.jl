# =================
# AlignmentAnchors
# =================


# Anchor definition
# ------------------
Docile.@doc """
A type to store the operation enocded in an alignment (CIGAR Operations).
Also stores the position in the alignment view of the sequences, and the
corresponding position in the unaltered source sequence (nucleotide or protein).

It stores the position values as Int types and the alignment operation is stored
as a type `Operation`, these are defined in the file `operations.jl` and loaded
into Bio.Align namespace as a series of global constants.

Calling the AlignmentAnchor default constructor with default arguments sets both
integer fields to 0, and the operation field to the global constant OP_INVALID.
""" ->
immutable AlignmentAnchor
  alnPos::Int
  srcPos::Int
  op::Operation
  function AlignmentAnchor(gp::Int = 0, sp::Int = 0, op::Operation = OP_INVALID)
    return new(gp, sp, op)
  end
end


# Basic operators for AlignmentAnchors
# -------------------------------------

Docile.@doc """
Print to screen or other IO stream, a formatted description of an
AlignmentAnchor.
""" ->
function show(io::IO, anc::AlignmentAnchor)
  write(io, "Alignment Position: $(anc.alnPos), Source Position: $(anc.srcPos), Operation: $(anc.op)")
end

Docile.@doc """
Copy an AlignmentAnchor to a new AlignmentAnchor variable.

Takes only one parameter of AlignmentAnchor type.
""" ->
function copy(src::AlignmentAnchor)
  return AlignmentAnchor(src.alnPos, src.srcPos, src.op)
end


# Basic operators for boolean operations consider
# positions, not operations

Docile.@doc """
Check for equity of two AlignmentAnchors.
""" ->
function ==(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.alnPos == b.alnPos && a.srcPos == b.srcPos
end

Docile.@doc """
Check for inequity of two AlignmentAnchors.
""" ->
function !=(a::AlignmentAnchor, b::AlignmentAnchor)
  return !(a == b)
end

Docile.@doc """
Check that AlignmentAnchor a, is less than AlignmentAnchor b.
""" ->
function <(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.alnPos < b.alnPos || a.srcPos < b.srcPos
end

Docile.@doc """
Check that AlignmentAnchor a, is greater than AlignmentAnchor b.
""" ->
function >(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.alnPos > b.alnPos || a.srcPos > b.srcPos
end

Docile.@doc """
Check that AlignmentAnchor a, is less than or equal to AlignmentAnchor b.
""" ->
function <=(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.alnPos <= b.alnPos || a.srcPos < b.srcPos
end

Docile.@doc """
Check that AlignmentAnchor a, is greater than or equal to AlignmentAnchor b.
""" ->
function >=(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.alnPos >= b.alnPos || a.srcPos > b.srcPos
end

Docile.@doc """
Check whether the alignment anchor contains the specified operation.
""" ->
function hasOp(anc::AlignmentAnchor, op::Operation)
    return anc.op == op
end



# AlignmentAnchors Definition
# ----------------------------

typealias AlignmentAnchors Vector{AlignmentAnchor}

function show(io::IO, aa::AlignmentAnchors)
  out::String = ""
  for i in aa
    out *= "($(i.srcPos), $(i.alnPos), $(Char(i.op)))"
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

Docile.@doc """
An immutable field-less type, inheriting from Ordering. This is used with
Julia's sorting API functions to determine how an array of AlignmentAnchors is
sorted.

Specifically, it specifies that the sorting API functions use less-than
methods functions that compare the srcPos fields of the AlignmentAnchors.
""" ->
immutable srcPosOrdering <: Ordering end

Docile.@doc """
An immutable field-less type, inheriting from Ordering. This is used with
Julia's sorting API functions to determine how an array of AlignmentAnchors is
sorted.

Specifically, it specifies that the sorting API functions use less-than
methods functions that compare the alnPos fields of the AlignmentAnchors.
""" ->
immutable alnPosOrdering <: Ordering end

const BY_SRC = srcPosOrdering()
const BY_ALN = alnPosOrdering()

# Is AlignmentAnchor a, less than AlignmentAnchor b, according to the source
# sequence co-ordinate.

Docile.@doc """
Returns true if the srcPos field of AlignmentAnchor a, is less than the srcPos
field of AlignmentAnchor b.
""" ->
function lt(o::srcPosOrdering, a::AlignmentAnchor, b::AlignmentAnchor)
  return a.srcPos < b.srcPos
end

Docile.@doc """
Returns true if the alnPos field of AlignmentAnchor a, is less than the alnPos
field of AlignmentAnchor b.
""" ->
function lt(o::alnPosOrdering, a::AlignmentAnchor, b::AlignmentAnchor)
  return a.alnPos < b.alnPos
end

Docile.@doc """
Returns true if the srcPos field of AlignmentAnchor a, is less than the Integer
b.
""" ->
function lt(o::srcPosOrdering, a::AlignmentAnchor, b::Int)
  return a.srcPos < b
end

Docile.@doc """
Returns true if the alnPos field of AlignmentAnchor a, is less than the Integer
b.
""" ->
function lt(o::alnPosOrdering, a::AlignmentAnchor, b::Int)
  return a.alnPos < b
end

Docile.@doc """
Returns true if the srcPos field of Integer a, is less than the srcPos field of
the AlignmentAnchor b.
""" ->
function lt(o::srcPosOrdering, a::Int, b::AlignmentAnchor)
  return a < b.srcPos
end

Docile.@doc """
Returns true if the srcPos field of Integer a, is less than the alnPos field of
the AlignmentAnchor b.
""" ->
function lt(o::alnPosOrdering, a::Int, b::AlignmentAnchor)
  return a < b.alnPos
end

function upperBoundAnchor(arr::AlignmentAnchors, i::Int, o::alnPosOrdering)
    return searchsortedlast(arr, i, o) + 1
end

function lowerBoundAnchor(arr::AlignmentAnchors, i::Int, o::alnPosOrdering)
    return searchsortedfirst(arr, i, o)
end


immutable AlignedSequence
    src # Sequence from Bio.seq.
    anchors::AlignmentAnchors
end


function findPosition(anchors::AlignmentAnchors, position::Int, ordering::alnPosOrdering)
    loBracket = searchsortedlast(anchors, position, ordering)
    hiBracket = loBracket + 1
    # We probably need code here to handle a few edge cases.
    return loBracket, hiBracket
end


function alnToSrc(alignedSeq::AlignedSequence, alnPosition::Int)
    lowIdx, highIdx = findPosition(alignedSeq, alnPosition, BY_ALN)
    highAnchor = alignedSeq.anchors[highIdx]
    lowAnchor = alignedSeq.anchors[lowIdx]
    sourceDistance = hiAnc.srcPos - loAnc.srcPos
    alignDistance = alnPosition - loAnc.alnPos
    if srcDist > alnDist
        sourcePosition = loAnc.srcPos + alnDist
    else
        sourcePosition = hiAnc.srcPos
    end
    return sourcePosition
end


function srcToAln(alignedSeq::AlignedSequence, alnPosition::Int)
    lowIdx, highIdx = findPosition(alignedSeq, alnPosition, BY_SRC)
end
