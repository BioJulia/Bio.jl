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
    stream::BufferedStreams.BufferedInputStream
    done::Bool
    header::SectionHeader
    record::Record
    nodes::Vector{BBI.RTreeLeafNode}
    current_node::Int
    n_records::UInt16
    current_record::UInt16
end

function Base.start(iter::OverlapIterator)
    nodes = BBI.find_overlapping_nodes(iter.reader.rtree, iter.chromid, iter.chromstart, iter.chromend)
    # dummy header
    stream = Libz.ZlibInflateInputStream(iter.reader.stream)
    header = SectionHeader(0, 0, 0, 0, 0, 0, 0, 0)
    return OverlapIteratorState(stream, false, header, Record(), nodes, 1, 0, 0)
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
        while state.current_record == state.n_records && state.current_node ≤ endof(state.nodes)
            node = state.nodes[state.current_node]
            seek(iter.reader.stream, node.data_offset)
            state.stream = Libz.ZlibInflateInputStream(iter.reader.stream)
            state.header = read(state.stream, SectionHeader)
            state.current_node += 1
            state.n_records = state.header.item_count
            state.current_record = 0
        end
        if state.current_record == state.n_records && state.current_node > endof(state.nodes)
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
    return record.header.chrom_id == chromid && !(record.chromend ≤ chromstart || record.chromstart ≥ chromend)
end
