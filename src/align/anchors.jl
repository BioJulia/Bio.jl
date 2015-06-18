# =================
# AlignmentAnchors
# =================

# Anchor definition
# ------------------

immutable AlignmentAnchor
  gapPos::Int
  seqPos::Int
  op::Operation
  function AlignmentAnchor(gp::Int = 0, sp::Int = 0, op::Operation = OP_INVALID)
    return new(gp, sp, op)
  end
end


# Basic operators for AlignmentAnchors
# -------------------------------------

function show(io::IO, anc::AlignmentAnchor)
  write(io, 
"AlignmentAnchor:
----------------
Gap-Space position: $(anc.gapPos)
Source-space position: $(anc.seqPos)
Alignment operation: $(anc.op)")
end


function copy(src::AlignmentAnchor)
  return AlignmentAnchor(src.gapPos, src.seqPos, src.op)
end


# Basic operators for boolean operations consider
# positions, not operations

function ==(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.gapPos == b.gapPos && a.seqPos == b.seqPos
end

function !=(a::AlignmentAnchor, b::AlignmentAnchor)
  return !(a == b)
end

function <(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.gapPos < b.gapPos || a.seqPos < b.seqPos
end

function >(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.gapPos > b.gapPos || a.seqPos > b.seqPos
end

function <=(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.gapPos <= b.gapPos || a.seqPos < b.seqPos
end

function >=(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.gapPos >= b.gapPos || a.seqPos > b.seqPos
end


# AlignmentAnchors Array
# -----------------------

typealias AlignmentAnchors Vector{AlignmentAnchor}

function show(io::IO, aa::AlignmentAnchors)
  out::String = ""
  for i in aa
    out *= "($(i.seqPos), $(i.gapPos), $(Char(i.o)))"
  end
  write(io, out)
end


#=
In order to make use of Julia's sorting API and avoid having to write our own
searchsortedfirst and searchsortedlast methods, whilst ALSO avoiding slowing down
things by passing functions about as arguments to other functions, we make use of
how the sort API constructs its searching and sorting algorithms. It uses types 
that inherit from the Ordering abstract type, and this type, combined with use of
different Base.Order.lt methods - define how the values in the array are compared
to the query value and in what order. So we take advantage of this and define our
own types inheriting from Base.Order.Ordering and create our own lt methods. These
will then be used by the core Julia sorting API efficiently.
=#

immutable seqPosOrdering <: Ordering end
immutable gapPosOrdering <: Ordering end


const BY_SEQ = seqPosOrdering()
const BY_GAP = gapPosOrdering()


function lt(o::seqPosOrdering, a::AlignmentAnchor, x::Int)
  return a.seqPos < x
end

function lt(o::gapPosOrdering, a::AlignmentAnchor, x::Int)
  return a.gapPos < x
end
