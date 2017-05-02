# BigBed Writer
# =============

type WriterState
    # section info
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32

    # last record info
    last_chrom_start::UInt32
    last_chrom_end::UInt32

    # section state
    started::Bool
    buffer::IOBuffer
    intervals::Vector{Interval{Void}}

    # global state
    finished_chrom_ids::Set{UInt32}
    summaries::Vector{BBI.SectionSummary}

    function WriterState()
        return new(
            # section info
            0, 0, 0,
            # last record info
            0, 0,
            # section state
            false, IOBuffer(), Interval{Void}[],
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

    # chrom ID => name
    chromnames::Dict{UInt32,String}

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
    chromlist_with_id = BBI.add_chrom_ids(chromlist)
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
        Dict(id => name for (name, id, len) in chromlist_with_id),
        WriterState())
end

#function Base.write(writer::Writer, record::Tuple{String,Integer,Integer,Vararg})
function Base.write(writer::Writer, record::Tuple{String,Integer,Integer})
    chromname, chromstart, chromend = record
    chromid, _ = writer.chroms[chromname]
    return write_impl(writer, chromid, UInt32(chromstart - 1), UInt32(chromend))
end

function Base.close(writer::Writer)
    state = writer.state
    if state.started
        finish_section!(writer)
    end

    # write data index
    stream = writer.stream
    data_index_offset = position(stream)
    BBI.write_rtree(writer.stream, state.summaries)

    # fill header
    header = BBI.Header(
        BBI.BED_MAGIC,
        3,  # version
        0,  # zoom level
        writer.chrom_tree_offset,
        writer.data_offset,
        data_index_offset,
        3,  # field count (0 for bigWig)
        3,  # predefined field count (0 for bigWig?)
        0,  # autoSQL offset (0 for bigWig?)
        writer.summary_offset,
        writer.uncompressed_buffer_size,
        0   # reserved
    )
    seekstart(stream)
    write(stream, header)

    # fill summary
    seek(stream, writer.summary_offset)
    write(stream, BBI.compute_total_sumamry(state.summaries))

    # fill data count
    seek(stream, writer.data_offset)
    write(stream, foldl(+, UInt64(0), s.count for s in state.summaries))

    close(stream)
    return
end

function write_impl(writer::Writer, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    # check consistency of new record
    if chromstart ≥ chromend
        throw(ArgumentError("non-positive interval"))
    end

    # TODO: fix
    state = writer.state
    if state.started && (chromid != state.chromid || position(state.buffer) ≥ writer.uncompressed_buffer_size - sizeof(UInt32) * 3)
        finish_section!(writer)
    end

    if !state.started
        start_section!(writer, chromid, chromstart, chromend)
    end

    # write data to the buffer
    n = write(state.buffer, chromid, chromstart, chromend, 0x00)
    state.chromend = max(state.chromend, chromend)

    # update last record info
    state.last_chrom_start = chromstart
    state.last_chrom_end = chromend

    push!(state.intervals, Interval(writer.chromnames[chromid], chromstart + 1, chromend))

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

    # initialize last record info
    state.last_chrom_start = 0
    state.last_chrom_end = 0

    # set flag
    state.started = true
    empty!(state.intervals)
    return
end

function finish_section!(writer::Writer)
    state = writer.state

    # write compressed data
    data = takebuf_array(state.buffer)
    offset = position(writer.stream)
    write(writer.stream, Libz.compress(data))

    # record section summary
    count, cov, minval, maxval, sumval, sumsqval = compute_section_summary(state.intervals)
    #@show count, cov, minval, maxval, sumval, sumsqval
    push!(
        state.summaries,
        BBI.SectionSummary(
            state.chromid,
            state.chromstart,
            state.chromend,
            offset,
            count,
            cov,
            minval,
            maxval,
            sumval,
            sumsqval))
    state.started = false
    return
end

function compute_section_summary(intervals::Vector{Interval{Void}})
    cov = 0
    minval = Inf
    maxval = -Inf
    sumval = sumsqval = 0.0
    for range in Bio.Intervals.coverage(intervals)
        depth = range.metadata
        size = range.last - range.first + 1
        cov += size
        minval = min(minval, depth)
        maxval = max(maxval, depth)
        sumval += depth * size
        sumsqval += depth^2 * size
    end
    return length(intervals), cov, minval, maxval, sumval, sumsqval
end

function write_zeros(stream::IO, n::Integer)
    for _ in 1:n
        write(stream, 0x00)
    end
    return n
end
