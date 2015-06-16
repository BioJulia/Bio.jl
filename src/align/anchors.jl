# =================
# AlignmentAnchors
# =================

# Anchor definition
# ------------------

immutable AlignmentAnchor
  gapPos::Int
  seqPos::Int
  op::Operation
  function AlignmentAnchor(gp = 0, sp = 0, op = OP_INVALID)
    return new(gp, sp, op)
  end
end


# Basic operators for AlignmentAnchors
# -------------------------------------

function copy(src::AlignmentAnchor)
  return GapAnchor(src.gapPos, src.seqPos, src.op)
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

# Sorting gap anchors by their sequence space position.
function sortSeqPos(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.seqPos < b.seqPos
end

# Sorting gap anchors by their gap space position.
function sortGapPos(a::AlignmentAnchor, b::AlignmentAnchor)
  return a.gapPos < b.gapPos
end


# AlignmentAnchors Array
# -----------------------

typealias AlignmentAnchors Vector{AlignmentAnchor}


function upperBound{T}(array::Vector{T}, first::Int, last::Int, value::T, comp = x -> value < x)
  count::Int = last - first
  idx::Int = 0
  step::Int = 0
  while count > 0
    idx = first
    step = Int(round(count / 2))
    idx += step
    if !comp(array[idx])
      idx += 1
      first = idx
      count -= step + 1
    else
      count = step
    end
  end
  return first
end
