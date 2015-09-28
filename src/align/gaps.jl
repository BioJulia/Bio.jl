abstract Gaps

# ArrayGaps are for more efficient operations with when row-wise access is more
# desirable and inserting / extending gaps is desired.

type ArrayGaps <: Gaps
    source
    # An array with alternating source/gap character counts. The array always starts with
    # source, so if the first 5 bases are a gap, then the array would start [0, 5, ...]
    array::Vector{Int}
    source_end_pos::Int
    source_unclipped_end_pos::Int
    clip::UnitRange{Int}

    function ArrayGaps(source, array = Int[], source_end_pos = 0, source_unclipped_end_pos = 0, clip_beign_pos = 0, clip_end_pos = 0)
        x = new()
        x.source = source
        x.array = array
        x.source_end_pos = source_end_pos
        x.source_unclipped_end_pos = source_unclipped_end_pos
        x.clip = clip_beign_pos:clip_end_pos
        return x
    end
end


function unclipped_length(x::ArrayGaps)
    return x.source_end_pos + x.source_unclipped_end_pos
end


function clear_gaps!(x::ArrayGaps)

end


function clear_clipping!(x::ArrayGaps)
    x.source_end_pos = length(x.source)
    x.source_unclipped_end_pos
    x.clip = 0:length(x.source)
end


function view_position(x::ArrayGaps, position::Int)
    viewposition = sourceposition = position - x.clip.start
    next_block = start(x.array)
    current_block_size, next_block = next(x.array, next_block)
    while true
        if done(x.array, next_block)
            return viewposition
        end
        viewposition += current_block_size
        current_block_size, next_block = next(x.array, next_block)
        if sourceposition < current_block_size
            return viewposition
        end
        sourceposition -= current_block_size
        current_block_size, next_block = next(x.array, next_block)
    end
end


function source_position(x::ArrayGaps, position::Int)
    # Initialise iteration through the data array.
    next_block = start(x.array)
    current_block_size, next_block = next(x.array, next_block)
    # Initialise variables representing the remaining distance from the source position.
    # And the distance covered across the source position.
    remaining = position
    covered = x.clip.start
    while true
        if done(x.array, next_block)
            return covered
        end
        if remaining <= current_block_size
            return covered
        end
        remaining -= current_block_size
        current_block_size, next_block = next(x.array, next_block)
        if remaining <= current_block_size
            return covered + remaining
        end
        covered += current_block_size
        remaining -= current_block_size
        current_block_size, next_block = next(x.array, next_block)
    end
end
