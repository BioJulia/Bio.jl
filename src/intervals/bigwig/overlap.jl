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

function eachoverlap(reader::Reader, interval::Bio.Intervals.Interval)
    if haskey(reader.chroms, interval.seqname)
        id, _ = reader.chroms[interval.seqname]
    else
        id = typemax(UInt32)
    end
    return OverlapIterator(reader, id, interval.first - 1, interval.last)
end

type OverlapIteratorState
    # inflating data stream
    stream::BufferedStreams.BufferedInputStream
    done::Bool
    header::SectionHeader
    record::Record
    offsets::Vector{UInt64}
    current_offset::Int
    n_records::UInt16
    current_record::UInt16
end

function Base.start(iter::OverlapIterator)
    offsets = BBI.find_overlapping_data_offsets(iter.reader.rtree, iter.chromid, iter.chromstart, iter.chromend)
    # dummy header and record
    stream = Libz.ZlibInflateInputStream(iter.reader.stream)
    header = SectionHeader(0, 0, 0, 0, 0, 0, 0, 0)
    record = Record()
    return OverlapIteratorState(stream, false, header, record, offsets, 1, 0, 0)
end

function Base.done(iter::OverlapIterator, state::OverlapIteratorState)
    advance!(iter, state)
    return state.done
end

function Base.next(iter::OverlapIterator, state::OverlapIteratorState)
    return state.record, state
end

function advance!(iter::OverlapIterator, state::OverlapIteratorState)
    while true
        while state.current_record == state.n_records && state.current_offset ≤ endof(state.offsets)
            offset = state.offsets[state.current_offset]
            seek(iter.reader.stream, offset)
            state.stream = Libz.ZlibInflateInputStream(iter.reader.stream)
            state.header = read(state.stream, SectionHeader)
            state.current_offset += 1
            state.n_records = state.header.item_count
            state.current_record = 0
        end
        if state.current_record == state.n_records && state.current_offset > endof(state.offsets)
            state.done = true
            return state
        end

        read_record!(state.stream, state)
        record = state.record
        if get(record.header).chrom_id == iter.chromid && !(record.chromend ≤ iter.chromstart || record.chromstart ≥ iter.chromend)
            return state
        end
    end
end

function read_record!(stream, state)
    header = state.header
    if isbedgraph(header)
        chromstart = read(stream, UInt32)
        chromend   = read(stream, UInt32)
    elseif isvarstep(header)
        chromstart = read(stream, UInt32)
        chromend   = chromstart + header.item_span
    elseif isfixedstep(header)
        chromstart = (state.current_record == 0 ? header.chrom_start : state.record.chromstart) + header.item_step
        chromend   = chromstart + header.item_span
    else
        throw(ArgumentError("invalid data type"))
    end
    value = read(stream, Float32)
    state.record = Record(header, chromstart, chromend, value)
    state.current_record += 1
    return state
end
