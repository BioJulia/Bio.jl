abstract Gaps

# ArrayGaps are for more efficient operations with when row-wise access is more
# desirable and inserting / extending gaps is desired.

type ArrayGaps <: Gaps
  source
  # An array with alternating source/gap character counts. The array always starts with
  # source, so if the first 5 bases are a gap, then the array would start [0, 5, ...]
  array::Vector{Int}
  sourceEndPos::Int
  sourceUnclippedEndPos::Int
  clip::UnitRange{Int}

  function ArrayGaps(source, array = Int[], sourceEndPos = 0, sourceUnclippedEndPos = 0, clipBeginPos = 0, clipEndPos = 0)
    x = new()
    x.source = source
    x.array = array
    x.sourceEndPos = sourceEndPos
    x.sourceUnclippedEndPos = sourceUnclippedEndPos
    x.clip = clipBeginPos:clipEndPos
    return x
  end
end



function unclippedLength(x::ArrayGaps)
  return x.sourceEndPos + x.sourceUnclippedEndPos
end

function clearGaps!(x::ArrayGaps)

end

function clearClipping!(x::ArrayGaps)
  x.sourceEndPos = length(x.source)
  x.sourceUnclippedEndPos
  x.clip = 0:length(x.source)
end

function sourcePosition(x::ArrayGaps, position::Int)
  # Initialise iteration through the data array.
  nextBlock = start(x.array)
  currentBlockSize, nextBlock = next(x.array, nextBlock)
  
  # Initialise variables representing the remaining distance from the source position.
  # And the distance covered across the source position. 
  remaining = position
  covered = x.clipBeginPos

  while true
    if done(x.array, nextBlock)
      return covered
    end
    if remaining <= currentBlockSize
      return sourcePos
    end
    remaining -= currentBlockSize
    currentBlockSize, nextBlock = next(x.array, nextBlock)
    if remaining <= currentBlockSize
      return sourcePos + viewPos
    end
    covered += currentBlockSize
    remaining -= currentBlockSize
    currentBlockSize, nextBlock = next(x.array, nextBlock)
  end

end

ACCTGAC
|||
ACC--TG-AC
1234445667 : SourcePos
1234567890 : ViewPos

pos = 6

viewpos = 6 , 3 , 1

sourcepos = 2
            


0,3,2,2,1,2
X X X

# AnchorGaps are more efficient for column wise / site wise access. A binary search can
# be done to find specific positions.

immutable GapAnchor
  gapPos::Int
  seqPos::Int
end

type AnchorGaps <: Gaps
  source
  anchors::Vector{GapAnchor}
  cutBegin::Int
  cutEnd::Int
  viewCutBegin::Int
  viewCutEnd::Int

  function AnchorGaps(source, anchors = GapAnchor[], cutBegin = 0, cutEnd = 0, viewCutBegin = 0, viewCutEnd = 0)
    x = new()
    x.source = source
    x.anchors = anchors
    x.cutBegin = cutBegin
    x.cutEnd
    x.viewCutBegin
    x.viewCutEnd
    return x
  end
end

function clearGaps(x::AnchorGaps)
  x.anchors = GapAnchor[]
  x.cutBegin = 0
  x.cutEnd = 0
  x.viewCutBegin = 0
  x.viewCutEnd = 0
end
