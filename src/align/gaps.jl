abstract Gaps

# ArrayGaps are for more efficient operations with when row-wise access is more
# desirable and inserting / extending gaps is desired.

type ArrayGaps{S <: Bio.Seq.Sequence} <: Gaps
    source::S
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

function viewPosition(x::ArrayGaps, position::Int)
    viewposition = sourceposition = position - x.clip.start
    nextBlock = start(x.array)
    currentBlockSize, nextBlock = next(x.array, nextBlock)
    while true
        if done(x.array, nextBlock)
            return viewposition
        end
        viewposition += currentBlockSize
        currentBlockSize, nextBlock = next(x.array, nextBlock)
        if sourceposition < currentBlockSize
            return viewposition
        end
        sourceposition -= currentBlockSize
        currentBlockSize, nextBlock = next(x.array, nextBlock)
    end
end


function sourcePosition(x::ArrayGaps, position::Int)
    # Initialise iteration through the data array.
    nextBlock = start(x.array)
    currentBlockSize, nextBlock = next(x.array, nextBlock)
    # Initialise variables representing the remaining distance from the source position.
    # And the distance covered across the source position.
    remaining = position
    covered = x.clip.start
    while true
        if done(x.array, nextBlock)
            return covered
        end
        if remaining <= currentBlockSize
            return covered
        end
        remaining -= currentBlockSize
        currentBlockSize, nextBlock = next(x.array, nextBlock)
        if remaining <= currentBlockSize
            return covered + remaining
        end
        covered += currentBlockSize
        remaining -= currentBlockSize
        currentBlockSize, nextBlock = next(x.array, nextBlock)
    end
end
