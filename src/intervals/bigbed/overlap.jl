# BigBed Overlap
# ==============

immutable OverlapIterator
    reader::Reader
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
end

function Base.eltype(::Type{OverlapIterator})
    return Record
end

function Base.iteratorsize(::Type{OverlapIterator})
    return Base.SizeUnknown()
end

function Bio.Intervals.eachoverlap(reader::Reader, interval::Bio.Intervals.Interval)
    if haskey(reader.chroms, interval.seqname)
        id, _ = reader.chroms[interval.seqname]
    else
        id = typemax(UInt32)
    end
    return OverlapIterator(reader, id, interval.first - 1, interval.last)
end

type OverlapIteratorState
    state::Bio.Ragel.State
    done::Bool
    record::Record
    offsets::Vector{UInt64}
    current_offset::Int
end

function Base.start(iter::OverlapIterator)
    offsets = BBI.find_overlapping_data_offsets(iter.reader.rtree, iter.chromid, iter.chromstart, iter.chromend)
    if !isempty(offsets)
        seek(iter.reader.stream, first(offsets))
    end
    return OverlapIteratorState(
        Bio.Ragel.State(
            data_machine.start_state,
            Libz.ZlibInflateInputStream(iter.reader.stream, reset_on_end=false)),
        isempty(offsets), Record(), offsets, isempty(offsets) ? 1 : 2)
end

function Base.done(iter::OverlapIterator, state::OverlapIteratorState)
    advance!(iter, state)
    return state.done
end

function Base.next(iter::OverlapIterator, state::OverlapIteratorState)
    return copy(state.record), state
end

function advance!(iter::OverlapIterator, state::OverlapIteratorState)
    while true
        while state.current_offset ≤ endof(state.offsets) && eof(state.state.stream)
            seek(iter.reader.stream, state.offsets[state.current_offset])
            state.state = Bio.Ragel.State(data_machine.start_state, Libz.ZlibInflateInputStream(iter.reader.stream, reset_on_end=false))
            state.current_offset += 1
        end
        if state.done || (state.current_offset > endof(state.offsets) && eof(state.state.stream))
            state.done = true
            return state
        end

        _read!(iter.reader, state.state, state.record)
        state.record.reader = iter.reader
        if overlaps(state.record, iter.chromid, iter.chromstart, iter.chromend)
            return state
        end
    end
end

function overlaps(record::Record, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    return record.chromid == chromid && !(record.chromend ≤ chromstart || record.chromstart ≥ chromend)
end
