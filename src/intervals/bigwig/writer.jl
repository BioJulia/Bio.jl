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
    count::UInt64

    # last record info
    last_chrom_start::UInt32
    last_chrom_end::UInt32

    # zoom info
    coverage::Vector{UInt32}
    minval::Vector{Float32}
    maxval::Vector{Float32}
    sumval::Vector{Float32}
    sumsqval::Vector{Float32}

    # section state
    started::Bool
    buffer::IOBuffer

    # global state
    finished_chrom_ids::Set{UInt32}
    summaries::Vector{BBI.SectionSummary}

    function WriterState(datatype::Symbol)
        return new(
            # section info
            encode_datatype(datatype), 0, 0, 0, 0, 0, 0,
            # last record info
            0, 0,
            # zoom info
            UInt32[], Float32[], Float32[], Float32[], Float32[],
            # section data
            false, IOBuffer(),
            # global info
            Set{UInt32}(), BBI.SectionSummary[])
    end
end

immutable Writer <: Bio.IO.AbstractWriter
    # output stream
    stream::IO

    # the size of the smallest zoom in base
    binsize::UInt32

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
end

"""
    BigWig.Writer(output::IO, chromlist)

Create a data writer of the bigWig file format.

Arguments
---------

* `output`: data sink
* `chromlist`: chromosome list with length

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
function Writer(output::IO, chromlist::Union{AbstractVector,Associative})
    # write dummy header (filled later)
    write_zeros(output, sizeof(BBI.Header))

    # TODO: write zoom header
    write_zeros(output, 0)

    # write dummy total summary (filled later)
    summary_offset = position(output)
    write_zeros(output, sizeof(BBI.Summary))

    # write chromosome B+-tree
    chrom_tree_offset = position(output)
    chromlist_with_id = BBI.add_chrom_ids(chromlist)
    BBI.write_btree(output, chromlist_with_id)

    # write dummy data count (filled later)
    data_offset = position(output)
    write_zeros(output, sizeof(UInt64))

    return Writer(
        output,
        128,  # binsize
        64 * 2^10,
        summary_offset,
        chrom_tree_offset,
        data_offset,
        Dict(name => (id, len) for (name, id, len) in chromlist_with_id),
        Dict(id => len for (name, id, len) in chromlist_with_id),
        WriterState(:bedgraph))
end

function Base.write(writer::Writer, record::Tuple{String,Integer,Integer,Real})
    chromname, chromstart, chromend, value = record
    chromid, _ = writer.chroms[chromname]
    return write_impl(writer, chromid, UInt32(chromstart - 1), UInt32(chromend), Float32(value))
end

function Base.close(writer::Writer)
    state = writer.state
    if state.started
        finish_section!(writer)
    end

    # write R-tree
    stream = writer.stream
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
    seek(stream, writer.summary_offset)
    write(stream, BBI.compute_total_sumamry(state.summaries))

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

    # udpate zoom info
    state.count += 1
    bin = div(chromstart, writer.binsize)
    binstart = bin * writer.binsize
    cov = BBI.coverage2((chromstart, chromend), (binstart, binstart + writer.binsize))
    state.coverage[bin+1] += cov
    state.minval[bin+1] = min(state.minval[bin+1], value)
    state.maxval[bin+1] = max(state.maxval[bin+1], value)
    # TODO: is this correct?
    meanval = value * writer.binsize / cov  # mean value per bin
    state.sumval[bin+1] += meanval
    state.sumsqval[bin+1] += meanval^2

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
    state.count = 0

    # initialize last record info
    state.last_chrom_start = 0
    state.last_chrom_end = 0

    # initialize zoom info
    bincount = cld(writer.chromlens[chromid], writer.binsize)
    state.coverage = zeros(UInt32, bincount)
    state.minval = Vector{Float32}(bincount)
    fill!(state.minval, Inf32)
    state.maxval = Vector{Float32}(bincount)
    fill!(state.maxval, -Inf32)
    state.sumval = zeros(Float32, bincount)
    state.sumsqval = zeros(Float32, bincount)

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
    data = Libz.compress(takebuf_array(state.buffer))
    offset = position(writer.stream)
    write(writer.stream, data)

    # record section summary
    push!(
        state.summaries,
        BBI.SectionSummary(
            state.chromid,
            state.chromstart,
            state.chromend,
            offset,
            sizeof(data)))
    state.started = false
    return
end

function write_zeros(stream::IO, n::Integer)
    for _ in 1:n
        write(stream, 0x00)
    end
    return n
end
