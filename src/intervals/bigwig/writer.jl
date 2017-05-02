# BigWig Writer
# =============

type WriterState
    # section info
    datatype::UInt8
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
    itemstep::UInt32
    itemspan::UInt32

    # last record info
    last_chrom_start::UInt32
    last_chrom_end::UInt32

    # section summary
    count::UInt64
    coverage::UInt32
    minval::Float32
    maxval::Float32
    sumval::Float32
    sumsqval::Float32

    # section data
    started::Bool
    buffer::IOBuffer

    # global info
    finished_chrom_ids::Set{UInt32}
    summaries::Vector{BBI.SectionSummary}

    function WriterState(datatype::Symbol)
        return new(
            # section info
            encode_datatype(datatype), 0, 0, 0, 0, 0,
            # last record info
            0, 0,
            # section summary
            0, 0, 0.0, 0.0, 0.0, 0.0,
            # section data
            false, IOBuffer(),
            # global info
            Set{UInt32}(), BBI.SectionSummary[])
    end
end

immutable Writer <: Bio.IO.AbstractWriter
    # output stream
    stream::IO

    # maximum size of uncompressed buffer
    uncompressed_buffer_size::UInt64

    # file offsets
    summary_offset::UInt64
    chrom_tree_offset::UInt64
    data_offset::UInt64

    # chrom name => (ID, length)
    chroms::Dict{String,Tuple{UInt32,UInt32}}

    # mutable state
    state::WriterState
end

function Writer(stream::IO, chromlist::Vector)
    # write dummy header (filled later)
    write_zeros(stream, sizeof(BBI.Header))

    # TODO: write zoom header
    write_zeros(stream, 0)

    # write dummy total summary (filled later)
    summary_offset = position(stream)
    write_zeros(stream, sizeof(BBI.Summary))

    # write chromosome B+-tree
    chrom_tree_offset = position(stream)
    chromlist_with_id = add_chrom_ids(chromlist)
    BBI.write_btree(stream, chromlist_with_id)

    # write dummy data count (filled later)
    data_offset = position(stream)
    write_zeros(stream, sizeof(UInt64))

    return Writer(
        stream,
        64 * 2^10,
        summary_offset,
        chrom_tree_offset,
        data_offset,
        Dict(name => (id, len) for (name, id, len) in chromlist_with_id),
        WriterState(:bedgraph))
end

function Base.write(writer::Writer, record::Tuple{String,Integer,Integer,Real})
    chromname, chromstart, chromend, value = record
    chromid, _ = writer.chroms[chromname]
    return write_impl(writer, chromid, UInt32(chromstart - 1), UInt32(chromend), Float32(value))
end

function Base.close(writer::Writer)
    state = writer.state
    stream = writer.stream
    if state.started > 0
        finish_section!(writer)
    end

    # write R-tree
    data_index_offset = position(stream)
    BBI.write_rtree(writer.stream, state.summaries)

    # fill header
    header = BBI.Header(
        BBI.WIG_MAGIC,
        3,  # version
        0,  # zoom level
        writer.chrom_tree_offset,
        writer.data_offset,
        data_index_offset,
        0,  # field count (0 for bigWig)
        0,  # predefined field count (0 for bigWig?)
        0,  # autoSQL offset (0 for bigWig?)
        writer.summary_offset,
        writer.uncompressed_buffer_size,
        0   # reserved
    )
    seekstart(stream)
    write(stream, header)

    # fill summary
    summary = BBI.compute_total_sumamry(state.summaries)
    seek(stream, writer.summary_offset)
    write(stream, summary)

    # fill data count
    seek(stream, writer.data_offset)
    write(stream, UInt64(length(state.summaries)))

    close(stream)
    return
end

function write_impl(writer::Writer, chromid::UInt32, chromstart::UInt32, chromend::UInt32, value::Float32)
    # check consistency of new record
    if chromstart ≥ chromend
        throw(ArgumentError("non-positive interval"))
    end

    state = writer.state
    if state.started && (chromid != state.chromid || position(state.buffer) ≥ writer.uncompressed_buffer_size - sizeof(UInt32) * 3)
        finish_section!(writer)
    end

    n = 0
    if !state.started
        n += start_section!(writer, chromid, chromstart, chromend)
    end

    # write data to the buffer
    if isbedgraph(state.datatype)
        n += write(state.buffer, chromstart, chromend, value)
    elseif isvarstep(state.datatype)
        n += write(state.buffer, chromstart, value)
    elseif isfixedstep(state.datatype)
        n += write(state.buffer, value)
    else
        assert(false)
    end
    state.chromend = max(state.chromend, chromend)

    # update last record info
    state.last_chrom_start = chromstart
    state.last_chrom_end = chromend

    # udpate section summary
    state.count += 1
    state.coverage += chromend - chromstart
    state.minval = min(state.minval, value)
    state.maxval = max(state.maxval, value)
    state.sumval += value
    state.sumsqval += value^2

    return n
end

function start_section!(writer::Writer, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    state = writer.state
    if state.started
        error("unfinished section")
    end

    # initialize section info
    state.chromid = chromid
    state.chromstart = chromstart
    state.chromend = chromend
    if isfixedstep(state.datatype) || isvarstep(state.datatype)
        # inter the span from the first record
        state.itemspan = chromend - chromstart
    end

    # initialize last record info
    state.last_chrom_start = 0
    state.last_chrom_end = 0

    # initialize section summary
    state.count = 0
    state.coverage = 0
    state.minval = Inf32
    state.maxval = -Inf32
    state.sumval = 0.0f0
    state.sumsqval = 0.0f0

    # write dummy section header (filled later)
    n = write_zeros(state.buffer, sizeof(SectionHeader))
    state.started = true
    return n
end

function finish_section!(writer::Writer)
    state = writer.state

    # fill section header
    seekstart(state.buffer)
    write(
        state.buffer,
        SectionHeader(
            state.chromid,
            state.chromstart,
            state.chromend,
            state.itemstep,
            state.itemspan,
            state.datatype,
            0x00,  # reserved
            state.count))

    # write compressed section
    data = takebuf_array(state.buffer)
    offset = position(writer.stream)
    write(writer.stream, Libz.compress(data))

    # record section summary
    push!(
        state.summaries,
        BBI.SectionSummary(
            state.chromid,
            state.chromstart,
            state.chromend,
            offset,
            state.count,
            state.coverage,
            state.minval,
            state.maxval,
            state.sumval,
            state.sumsqval))
    state.started = false
    return
end

# Add a unique ID for each chromosome; chromlist is a tuple of (chrom name, chrom length).
function add_chrom_ids{T<:Integer}(chromlist::Vector{Tuple{String,T}})
    return [(name, UInt32(id - 1), UInt32(len)) for (id, (name, len)) in enumerate(chromlist)]
end

function write_zeros(stream::IO, n::Integer)
    for _ in 1:n
        write(stream, 0x00)
    end
    return n
end
