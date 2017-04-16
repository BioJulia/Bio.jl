# BigWig Reader
# =============

immutable Reader <: Bio.IO.AbstractReader
    stream::BufferedStreams.BufferedInputStream
    header::BBI.Header
    summary::BBI.TotalSummary
    btree::BBI.BTree
    rtree::BBI.RTree
    # chrom name => (ID, length)
    chroms::Dict{String,Tuple{UInt32,Int}}
end

function Base.eltype(::Type{Reader})
    return Record
end

function Reader(stream::IO)
    return Reader(BufferedStreams.BufferedInputStream(stream))
end

function Reader(stream::BufferedStreams.BufferedInputStream)
    # read header
    header = read(stream, BBI.Header)
    if header.magic != BBI.WIG_MAGIC
        error("invalid bigWig magic")
    elseif header.version < 3
        error("not a supported version of bigWig")
    end
    # read summary, B tree, and R tree
    seek(stream, header.total_summary_offset)
    summary = read(stream, BBI.TotalSummary)
    btree = BBI.BTree(stream, header.chromosome_tree_offset)
    rtree = BBI.RTree(stream, header.full_index_offset)
    chroms = Dict(name => (id, Int(len)) for (name, id, len) in BBI.chromlist(btree))
    return Reader(stream, header, summary, btree, rtree, chroms)
end


# Record
# ------

type Record
    header::Nullable{SectionHeader}
    chromstart::UInt32
    chromend::UInt32
    value::Float32
end


# Iterator
# --------

type IteratorState
    stream::BufferedStreams.BufferedInputStream
    done::Bool
    header::SectionHeader
    record::Record
    n_sections::UInt64
    current_section::UInt64
    n_records::UInt16
    current_record::UInt16
end

function Base.start(reader::Reader)
    seek(reader.stream, reader.header.full_data_offset)
    # this is defined as UInt32 in the spces but actually UInt64
    section_count = read(reader.stream, UInt64)
    # dummy header and record
    header = SectionHeader(0, 0, 0, 0, 0, 0, 0, 0)
    record = Record()
    state = IteratorState(Libz.ZlibInflateInputStream(reader.stream), false, header, record, section_count, 0, 0, 0)
    advance!(state)
    return state
end

function Base.done(reader::Reader, state)
    advance!(state)
    return state.done
end

function Base.next(reader::Reader, state)
    return state.record, state
end

function advance!(state::IteratorState)
    # find a section that has at least one record
    while state.current_section < state.n_sections && state.current_record == state.n_records
        state.header = read(state.stream, SectionHeader)
        state.current_section += 1
        state.n_records = state.header.item_count
        state.current_record = 0
    end
    if state.current_record == state.n_records
        state.done = true
        return state
    end
    @assert !state.done

    # read a new record
    header = state.header
    if isbedgraph(header)
        chromstart = read(state.stream, UInt32)
        chromend   = read(state.stream, UInt32)
    elseif isvarstep(header)
        chromstart = read(state.stream, UInt32)
        chromend   = chromstart + header.item_span
    elseif isfixedstep(header)
        chromstart = (state.current_record == 0 ? header.chrom_start : state.record.chromstart) + header.item_step
        chromend   = chromstart + header.item_span
    else
        throw(ArgumentError("invalid data type"))
    end
    value = read(state.stream, Float32)
    state.record = Record(header, chromstart, chromend, value)
    state.current_record += 1
    return state
end
