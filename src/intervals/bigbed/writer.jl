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
    compressed::Vector{UInt8}
    intervals::Vector{Interval{Void}}

    # record buffer
    recordbuffer::IOBuffer

    # global state and stats
    blocks::Vector{BBI.Block}
    finished_chrom_ids::Set{UInt32}
    nfields::Int
    count::UInt64
    cov::UInt32
    min::Float32
    max::Float32
    sum::Float32
    ssq::Float32

    function WriterState()
        return new(
            # section info
            0, 0, 0,
            # last record info
            0, 0,
            # section state
            false, IOBuffer(), UInt8[], Interval{Void}[],
            # record buffer
            IOBuffer(),
            # global info
            BBI.Block[], Set{UInt32}(), 0, 0,
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

    # chrom ID => name
    chromnames::Dict{UInt32,String}

    # mutable state
    state::WriterState

    # zoom buffer
    zoombuffer::BBI.ZoomBuffer
end

const ZOOM_SCALE = 4

"""
    BigBed.Writer(output::IO, chromlist; binsize=64)

Create a data writer of the bigBed file format.

Arguments
---------

* `output`: data sink
* `chromlist`: chromosome list with length
* `binsize=64`: size of a zoom with the highest resolution

Examples
--------

```julia
output = open("data.bb", "w")
writer = BigBed.Writer(output, [("chr1", 12345), ("chr2", 9100)])
write(writer, ("chr1", 101, 150, "gene 1"))
write(writer, ("chr2", 211, 250, "gene 2"))
close(writer)
```
"""
function Writer(output::IO, chromlist::Union{AbstractVector,Associative};
                binsize::Integer=64)
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
    max_block_size = 64 * 2^10
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
        Dict(id => name for (name, id, len) in chromlist_with_id),
        WriterState(),
        zoombuffer)
end

function Base.write(writer::Writer, record::Tuple{String,Integer,Integer,Vararg})
    chromname, chromstart, chromend = record
    optionals = Base.tail(Base.tail(Base.tail(record)))
    chromid, _ = writer.chroms[chromname]
    return write_impl(writer, chromid, UInt32(chromstart - 1), UInt32(chromend), optionals)
end

function Base.close(writer::Writer)
    state = writer.state
    if state.started
        finish_section!(writer)
    end

    # write data index
    stream = writer.stream
    data_index_offset = position(stream)
    BBI.write_rtree(writer.stream, state.blocks)

    # fill header
    header = BBI.Header(
        BBI.BED_MAGIC,
        3,  # version
        writer.zoomlevels,
        writer.chrom_tree_offset,
        writer.data_offset,
        data_index_offset,
        state.nfields,
        state.nfields,  # predefined field count (should be same as above?)
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
    write(stream, state.count)

    # write zoom
    seekend(stream)
    zoomheaders = BBI.write_zoom(stream, writer.zoombuffer, writer.zoomlevels, ZOOM_SCALE)

    # fill zoom headers
    seek(stream, BBI.HEADER_SIZE)
    for zheader in zoomheaders
        write(stream, zheader)
    end

    close(stream)
    return
end

function write_impl(writer::Writer, chromid::UInt32, chromstart::UInt32, chromend::UInt32, optionals::Tuple)
    # check consistency of new record
    state = writer.state
    if chromstart ≥ chromend
        throw(ArgumentError("non-positive interval"))
    end
    if state.nfields == 0
        # infer the number of fields from the first record
        state.nfields = 3 + length(optionals)
    elseif state.nfields != 3 + length(optionals)
        throw(ArgumentError("inconsistent field counts"))
    end

    # write record to record buffer in order to estimate the section buffer size
    truncate(state.recordbuffer, 0)
    write(state.recordbuffer, chromid, chromstart, chromend)
    write_optionals(state.recordbuffer, optionals)
    write(state.recordbuffer, 0x00)
    if state.started && (chromid != state.chromid || position(state.buffer) + position(state.recordbuffer) > writer.uncompressed_buffer_size)
        finish_section!(writer)
    end

    if !state.started
        start_section!(writer, chromid, chromstart, chromend)
    end

    # write data to the buffer
    seekstart(state.recordbuffer)
    n = write(state.buffer, state.recordbuffer)
    state.chromend = max(state.chromend, chromend)

    # update last record info
    state.last_chrom_start = chromstart
    state.last_chrom_end = chromend

    push!(state.intervals, Interval(writer.chromnames[chromid], chromstart + 1, chromend))

    return n
end

# Write optional fields to stream.
function write_optionals(stream::IO, optionals::Tuple)
    n = length(optionals)
    if n ≥ 1  # name
        print(stream, optionals[1])
    end
    if n ≥ 2  # score
        print(stream, '\t', optionals[2])
    end
    if n ≥ 3  # strand
        print(stream, '\t', optionals[3])
    end
    if n ≥ 4  # thickstart
        print(stream, '\t', optionals[4] - 1)
    end
    if n ≥ 5  # thickend
        print(stream, '\t', optionals[5])
    end
    if n ≥ 6  # itemrgb
        print(stream, '\t')
        print_rgb(stream, optionals[6])
    end
    if n ≥ 7  # blockcount
        print(stream, '\t', optionals[7])
    end
    if n ≥ 8  # blocksizes
        print(stream, '\t')
        for x in optionals[8]
            print(stream, x, ',')
        end
    end
    if n ≥ 9  # blockstarts
        print(stream, '\t')
        for x in optionals[9]
            print(stream, x - 1, ',')
        end
    end
end

function print_rgb(stream::IO, rgb::ColorTypes.RGB)
    print_rgb(stream, convert(ColorTypes.RGB{N0f8}, rgb))
end

function print_rgb(stream::IO, rgb::ColorTypes.RGB{N0f8})
    print(stream, reinterpret(UInt8, rgb.r), ',', reinterpret(UInt8, rgb.g), ',', reinterpret(UInt8, rgb.b))
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

    # initialize data buffer
    truncate(state.buffer, 0)
    resize!(state.compressed, div(writer.uncompressed_buffer_size * 11, 10))

    # set flag
    state.started = true
    empty!(state.intervals)
    return
end

function finish_section!(writer::Writer)
    state = writer.state

    # write compressed data
    datasize = BBI.compress!(state.compressed, takebuf_array(state.buffer))
    offset = position(writer.stream)
    unsafe_write(writer.stream, pointer(state.compressed), datasize)

    # record block
    push!(
        state.blocks,
        BBI.Block((state.chromid, state.chromstart), (state.chromid, state.chromend), offset, datasize))

    # update global stats
    count, cov, min, max, sum, ssq = compute_section_summary(state.intervals)
    state.count += count
    state.cov += cov
    state.min = Base.min(state.min, min)
    state.max = Base.max(state.max, max)
    state.sum += sum
    state.ssq += ssq

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
