# BigWig Writer
# =============

type WriterState
    # section info
    datatype::UInt8
    count::UInt64
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
    itemstep::UInt32
    itemspan::UInt32

    # previous record info
    chromstart_prev::UInt32
    chromend_prev::UInt32

    # section state
    started::Bool
    buffer::IOBuffer
    compressed::Vector{UInt8}

    # global state and stats
    blocks::Vector{BBI.Block}
    finished_chrom_ids::Set{UInt32}
    cov::UInt32
    min::Float32
    max::Float32
    sum::Float32
    ssq::Float32

    function WriterState(datatype::Symbol)
        return new(
            # section info
            encode_datatype(datatype), 0, 0, 0, 0, 0, 0,
            # last record info
            0, 0,
            # section data
            false, IOBuffer(), UInt8[],
            # global info
            BBI.Block[], Set{UInt32}(),
            0, +Inf32, -Inf32, 0.0f0, 0.0f0)
    end
end

immutable Writer <: Bio.IO.AbstractWriter
    # output stream
    stream::IO

    # number of zoom levels
    zoomlevels::Int

    # maximum size of uncompressed buffer
    uncompressed_buffer_size::UInt64

    # file offsets
    summary_offset::UInt64
    chrom_tree_offset::UInt64
    data_offset::UInt64

    # chrom name => (ID, length)
    chroms::Dict{String,Tuple{UInt32,UInt32}}

    # chrom ID => length
    chromlens::Dict{UInt32,UInt32}

    # mutable state
    state::WriterState

    # zoom buffer
    zoombuffer::BBI.ZoomBuffer
end

const ZOOM_SCALE = 4

"""
    BigWig.Writer(output::IO, chromlist; binsize=64, datatype=:bedgraph)

Create a data writer of the bigWig file format.

Arguments
---------

* `output`: data sink
* `chromlist`: chromosome list with length
* `binsize=64`: size of a zoom with the highest resolution
* `datatype=:bedgraph`: encoding of values (`:bedgraph`, `:varstep` or `:fixedstep`)

Examples
--------

```julia
output = open("data.bw", "w")
writer = BigWig.Writer(output, [("chr1", 12345), ("chr2", 9100)])
write(writer, ("chr1", 501, 600, 1.0))
write(writer, ("chr2", 301, 450, 3.0))
close(writer)
```
"""
function Writer(output::IO, chromlist::Union{AbstractVector,Associative};
                binsize::Integer=64, datatype::Symbol=:bedgraph)
    # write dummy header (filled later)
    write_zeros(output, BBI.HEADER_SIZE)

    # write dummy zoom headers (filled later)
    chromlist_with_id = BBI.add_chrom_ids(chromlist)
    maxlen = Base.maximum(x[3] for x in chromlist_with_id)
    zoomlevels = BBI.determine_zoomlevels(maxlen, binsize, ZOOM_SCALE)
    write_zeros(output, BBI.ZOOM_HEADER_SIZE * zoomlevels)

    # write dummy total summary (filled later)
    summary_offset = position(output)
    write_zeros(output, BBI.SUMMARY_SIZE)

    # write chromosome B+-tree
    chrom_tree_offset = position(output)
    BBI.write_btree(output, chromlist_with_id)

    # write dummy data count (filled later)
    data_offset = position(output)
    write_zeros(output, sizeof(UInt64))

    # initialize zoom buffer (use temporary file?)
    max_block_size = 16 * 2^10
    chromlens = Dict(id => len for (name, id, len) in chromlist_with_id)
    zoombuffer = BBI.ZoomBuffer(chromlens, binsize, max_block_size)

    return Writer(
        output,
        zoomlevels,
        max_block_size,
        summary_offset,
        chrom_tree_offset,
        data_offset,
        Dict(name => (id, len) for (name, id, len) in chromlist_with_id),
        chromlens,
        WriterState(datatype),
        zoombuffer)
end

function Base.write(writer::Writer, record::Tuple{String,Integer,Integer,Real})
    chromname, chromstart, chromend, value = record
    chromid = writer.chroms[chromname][1]
    return write_impl(writer, chromid, UInt32(chromstart - 1), UInt32(chromend), Float32(value))
end

function Base.write{T<:Real}(writer::Writer, interval::Bio.Intervals.Interval{T})
    chromid = writer.chroms[interval.seqname][1]
    return write_impl(writer, chromid, UInt32(interval.first - 1), UInt32(interval.last), Float32(interval.metadata))
end

function Base.close(writer::Writer)
    state = writer.state
    if state.started
        finish_section!(writer)
    end

    # write R-tree
    stream = writer.stream
    data_index_offset = position(stream)
    BBI.write_rtree(writer.stream, state.blocks)

    # fill header
    header = BBI.Header(
        BBI.WIG_MAGIC,
        3,  # version
        writer.zoomlevels,
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
    seek(stream, writer.summary_offset)
    write(stream, BBI.Summary(state.cov, state.min, state.max, state.sum, state.ssq))

    # fill data count
    seek(stream, writer.data_offset)
    write(stream, UInt64(length(state.blocks)))

    # write zoom
    seekend(stream)
    zoomheaders = BBI.write_zoom(stream, writer.zoombuffer, writer.zoomlevels, ZOOM_SCALE)
    @assert length(zoomheaders) == writer.zoomlevels

    # fill zoom headers
    seek(stream, BBI.HEADER_SIZE)
    for zheader in zoomheaders
        write(stream, zheader)
    end

    close(stream)
    return
end

function Bio.IO.stream(writer::Writer)
    return writer.stream
end

function write_impl(writer::Writer, chromid::UInt32, chromstart::UInt32, chromend::UInt32, value::Float32)
    n = 0
    state = writer.state
    if state.started && (chromid != state.chromid || position(state.buffer) ≥ writer.uncompressed_buffer_size - sizeof(UInt32) * 3)
        finish_section!(writer)
    end
    if !state.started
        n += start_section!(writer, chromid, chromstart, chromend)
    end

    # check new record
    check_interval(state, chromid, chromstart, chromend)

    # infer the step size from the first and second records
    if isfixedstep(state.datatype) && state.count == 1
        state.itemstep = chromstart - state.chromstart_prev
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

    # update state
    state.count += 1
    state.chromend = max(state.chromend, chromend)
    state.chromstart_prev = chromstart
    state.chromend_prev = chromend

    size = chromend - chromstart
    state.cov += size
    state.min = min(state.min, value)
    state.max = max(state.max, value)
    state.sum += value   * size
    state.ssq += value^2 * size

    BBI.add_value!(writer.zoombuffer, chromid, chromstart, chromend, value)

    return n
end

function check_interval(state::WriterState, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    if chromstart ≥ chromend
        throw(ArgumentError("empty interval"))
    elseif BBI.compare_intervals((chromid, chromstart, chromend), (state.chromid, state.chromstart_prev, state.chromend_prev)) != 1
        throw(ArgumentError("disordered intervals"))
    elseif isvarstep(state.datatype) && chromend - chromstart != state.itemspan
        throw(ArgumentError("inconsistent interval span"))
    elseif isfixedstep(state.datatype)
        if chromend - chromstart != state.itemspan
            throw(ArgumentError("inconsistent interval span"))
        elseif state.count > 1 && chromstart - state.chromstart_prev != state.itemstep
            throw(ArgumentError("inconsistent intreval step"))
        end
    end
end

function start_section!(writer::Writer, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    state = writer.state
    @assert !state.started

    # initialize section info
    state.chromid = chromid
    state.chromstart = chromstart
    state.chromend = chromend
    if isfixedstep(state.datatype) || isvarstep(state.datatype)
        # inter the span from the first record
        state.itemspan = chromend - chromstart
    end
    state.count = 0

    # initialize last record info
    state.chromstart_prev = 0
    state.chromend_prev = 0

    # write dummy section header (filled later)
    seekstart(state.buffer)
    truncate(state.buffer, SECTION_HEADER_SIZE)
    seek(state.buffer, SECTION_HEADER_SIZE)
    resize!(state.compressed, div(writer.uncompressed_buffer_size * 11, 10))
    state.started = true
    return SECTION_HEADER_SIZE
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
    size = BBI.compress!(state.compressed, takebuf_array(state.buffer))
    offset = position(writer.stream)
    unsafe_write(writer.stream, pointer(state.compressed), size)

    # record block
    push!(
        state.blocks,
        BBI.Block((state.chromid, state.chromstart), (state.chromid, state.chromend), offset, size))
    state.started = false
    return
end

function write_zeros(stream::IO, n::Integer)
    for _ in 1:n
        write(stream, 0x00)
    end
    return n
end
