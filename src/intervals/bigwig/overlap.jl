# BigWig Overlap
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
    # inflating data stream
    stream::IOBuffer
    data::Vector{UInt8}
    done::Bool
    header::SectionHeader
    record::Record
    blocks::Vector{BBI.Block}
    current_block::Int
    n_records::UInt16
    current_record::UInt16
end

function Base.start(iter::OverlapIterator)
    data = Vector{UInt8}(iter.reader.header.uncompress_buf_size)
    blocks = BBI.find_overlapping_blocks(iter.reader.index, iter.chromid, iter.chromstart, iter.chromend)
    # dummy header
    header = SectionHeader(0, 0, 0, 0, 0, 0, 0, 0)
    return OverlapIteratorState(IOBuffer(), data, false, header, Record(), blocks, 1, 0, 0)
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
        # find a section that has at least one record
        while state.current_record == state.n_records && state.current_block ≤ endof(state.blocks)
            block = state.blocks[state.current_block]
            seek(iter.reader.stream, block.offset)
            size = BBI.uncompress!(state.data, read(iter.reader.stream, block.size))
            state.stream = IOBuffer(state.data[1:size])
            state.header = read(state.stream, SectionHeader)
            state.current_block += 1
            state.n_records = state.header.itemcount
            state.current_record = 0
        end
        if state.current_record == state.n_records && state.current_block > endof(state.blocks)
            state.done = true
            return state
        end

        # read a new record
        _read!(iter.reader, state, state.record)
        if overlaps(state.record, iter.chromid, iter.chromstart, iter.chromend)
            return state
        end
    end
end

function overlaps(record::Record, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    return record.header.chromid == chromid && !(record.chromend ≤ chromstart || record.chromstart ≥ chromend)
end
