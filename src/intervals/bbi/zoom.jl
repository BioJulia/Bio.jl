# BBI Zoom Header/Data
# ====================

# Supplemental Table 6.
immutable ZoomHeader
    reduction_level::UInt32
    reserved::UInt32
    dataoffset::UInt64
    indexoffset::UInt64
end

const ZOOM_HEADER_SIZE = 24

function Base.read(stream::IO, ::Type{ZoomHeader})
    u32() = read(stream, UInt32)
    u64() = read(stream, UInt64)
    return ZoomHeader(u32(), u32(), u64(), u64())
end

function Base.write(stream::IO, header::ZoomHeader)
    return write(
        stream,
        header.reduction_level,
        header.reserved,
        header.dataoffset,
        header.indexoffset)
end

immutable Zoom{T<:IO}
    header::ZoomHeader
    rtree::RTree{T}
    # preallocated ZoomData buffer
    data::Vector{UInt8}
end

function Zoom(stream::IO, header::ZoomHeader, maxsize::Integer)
    return Zoom(header, RTree(stream, header.indexoffset), Vector{UInt8}(maxsize))
end

# Supplemental Table 19.
immutable ZoomData
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
    cov::UInt32
    min::Float32
    max::Float32
    sum::Float32
    ssq::Float32
end

const ZOOM_DATA_SIZE = 32

function Base.read(stream::IO, ::Type{ZoomData})
    u32() = read(stream, UInt32)
    f32() = read(stream, Float32)
    return ZoomData(
        u32(), u32(), u32(), u32(),
        f32(), f32(), f32(), f32())
end

function Base.write(stream::IO, data::ZoomData)
    return write(
        stream,
        data.chromid, data.chromstart, data.chromend,
        data.cov, data.min, data.max, data.sum, data.ssq)
end

function find_overlapping_zoomdata(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    blocks = find_overlapping_blocks(zoom.rtree, chromid, chromstart, chromend)
    stream = zoom.rtree.stream
    ret = ZoomData[]
    for block in blocks
        seek(stream, block.offset)
        size = uncompress!(zoom.data, read(stream, block.size))
        datastream = IOBuffer(zoom.data[1:size])
        while !eof(datastream)
            data = read(datastream, ZoomData)
            c = compare_intervals((data.chromid, data.chromstart, data.chromend), (chromid, chromstart, chromend))
            if c == 0  # overlap
                push!(ret, data)
            elseif c > 0
                break
            end
        end
    end
    return ret
end

function compare_intervals(x::NTuple{3,UInt32}, y::NTuple{3,UInt32})
    if x[1] < y[1] || (x[1] == y[1] && x[3] ≤ y[2])
        # strictly left
        return -1
    elseif x[1] > y[1] || (x[1] == y[1] && x[2] ≥ y[3])
        # strictly right
        return +1
    else
        # overlapping
        return 0
    end
end


# Statistics
# ----------

function coverage(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    count::UInt32 = 0
    if chromstart ≥ chromend  # empty range
        return count
    end
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    for d in data
        cov = coverage2((d.chromstart, d.chromend), (chromstart, chromend))
        count += round(UInt32, d.cov * cov / (d.chromend - d.chromstart))
    end
    return count
end

function mean(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    sum = 0.0f0
    size = 0
    for d in data
        cov = coverage2((d.chromstart, d.chromend), (chromstart, chromend))
        sum += (d.sum * cov) / (d.chromend - d.chromstart)
        size += cov
    end
    return sum / size
end

function minimum(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    return isempty(data) ? NaN32 : foldl((x, d) -> min(x, d.min), +Inf32, data)
end

function maximum(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    return isempty(data) ? NaN32 : foldl((x, d) -> max(x, d.max), -Inf32, data)
end

function std(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    sum = 0.0f0
    ssq = 0.0f0
    size = 0
    for d in data
        cov = coverage2((d.chromstart, d.chromend), (chromstart, chromend))
        sum += (d.sum * cov) / (d.chromend - d.chromstart)
        ssq += (d.ssq * cov) / (d.chromend - d.chromstart)
        size += cov
    end
    return sqrt((ssq - sum^2 / size) / (size - 1))
end

# Find the best zoom level for a given size.
function find_best_zoom(zooms::Vector{Zoom}, size::UInt32)::Nullable{Zoom}
    # NOTE: This assumes zooms are sorted by reduction_level.
    halfsize = div(size, 2)
    i = 0
    while i ≤ endof(zooms) && zooms[i+1].header.reduction_level ≤ halfsize
        i += 1
    end
    if i == 0
        return Nullable{Zoom}()
    else
        return Nullable(zooms[i])
    end
end

# Compute the nubmer of overlapping bases of [x1, x2) and [y1, y2).
function coverage2(x::NTuple{2,UInt32}, y::NTuple{2,UInt32})
    x1, x2 = x
    y1, y2 = y
    return max(UInt32(0), min(x2, y2) - max(x1, y1))
end

# Determine the minimum number of zoom levels that covers `len`.
function determine_zoomlevels(len::Integer, binsize::Integer, scale::Integer)
    if binsize < 0 || scale < 1
        return 0
    end
    level = 1
    while binsize * scale^level < len
        level += 1
    end
    return level
end


# Zoom Buffer
# -----------

type ZoomBuffer
    # chrom ID => length
    chromlens::Dict{UInt32,UInt32}

    # bin size of zoom
    binsize::UInt32

    # maximum size of uncompressed buffer
    max_buffer_size::UInt64

    # current chromosome (typemax(UInt32) for undefined)
    chromid::UInt32

    # zoom blocks
    blocks::Vector{Block}

    # zoom statistics
    cov::Vector{UInt32}
    min::Vector{Float32}
    max::Vector{Float32}
    sum::Vector{Float32}
    ssq::Vector{Float32}

    # flag whether the zoom statistics are buffered
    buffered::Bool

    # temporary files for data and directory
    stream1::IOStream
    stream2::IOStream
    tmpdir::String

    # number of zoom records written in stream1
    count::Int

    function ZoomBuffer(chromlens, binsize, max_block_size)
        tmpdir = mktempdir()
        try
            stream1 = mktemp(tmpdir)[2]
            stream2 = mktemp(tmpdir)[2]
            buffer = new(
                chromlens,
                binsize,
                max_block_size,
                typemax(UInt32),
                Block[],
                UInt32[], Float32[], Float32[], Float32[], Float32[],
                false,
                stream1, stream2, tmpdir,
                0)
            finalizer(buffer, buffer -> rm(buffer.tmpdir; force=true, recursive=true))
            return buffer
        catch
            rm(tmpdir; recursive=true)
            rethrow()
        end
    end
end

function add_value!(buffer::ZoomBuffer, chromid::UInt32, chromstart::UInt32, chromend::UInt32, value::Float32)
    if buffer.chromid == typemax(UInt32)
        init_buffer!(buffer, chromid)
    elseif chromid != buffer.chromid
        write_buffered_data(buffer)
        init_buffer!(buffer, chromid)
    end
    @assert chromid == buffer.chromid
    binsize = buffer.binsize
    for bin in div(chromstart, binsize):div(chromend - 1, binsize)
        binstart = UInt32(bin * binsize)
        cov = coverage2((binstart, binstart + binsize), (chromstart, chromend))
        bin1 = bin + 1
        buffer.cov[bin1] += cov
        buffer.min[bin1]  = min(buffer.min[bin1], value)
        buffer.max[bin1]  = max(buffer.max[bin1], value)
        buffer.sum[bin1] += value   * cov
        buffer.ssq[bin1] += value^2 * cov
    end
    buffer.buffered = true
    return buffer
end

function init_buffer!(buffer::ZoomBuffer, chromid::UInt32)
    @assert !buffer.buffered
    nbins = cld(buffer.chromlens[chromid], buffer.binsize)
    resize!(buffer.cov, nbins)
    fill!(buffer.cov, 0)
    resize!(buffer.min, nbins)
    fill!(buffer.min, +Inf32)
    resize!(buffer.max, nbins)
    fill!(buffer.max, -Inf32)
    resize!(buffer.sum, nbins)
    fill!(buffer.sum, 0.0f0)
    resize!(buffer.ssq, nbins)
    fill!(buffer.ssq, 0.0f0)
    buffer.chromid = chromid
    return buffer
end

function write_buffered_data(buffer::ZoomBuffer)
    @assert buffer.chromid != typemax(UInt32)
    # write zoom data to a temporary file
    for i in 1:endof(buffer.cov)
        chromstart = (i - 1) * buffer.binsize
        write(
            buffer.stream1,
            ZoomData(
                buffer.chromid,
                chromstart,
                chromstart + buffer.binsize,
                buffer.cov[i],
                buffer.min[i],
                buffer.max[i],
                buffer.sum[i],
                buffer.ssq[i]))
        buffer.count += 1
    end
    buffer.buffered = false
    return
end

function write_zoom(output::IO, buffer::ZoomBuffer, nlevels::Int, scale::Int)
    if buffer.buffered
        write_buffered_data(buffer)
    end
    headers = ZoomHeader[]
    for l in 1:nlevels
        binsize = buffer.binsize * scale^(l - 1)
        dataoffset, indexoffset = write_zoom_impl(output, buffer, scale)
        push!(headers, ZoomHeader(binsize, 0, dataoffset, indexoffset))
    end
    return headers
end

function write_zoom_impl(output::IO, buffer::ZoomBuffer, scale::Int)
    lo = (typemax(UInt32), typemax(UInt32))
    up = (typemin(UInt32), typemin(UInt32))
    blocks = Block[]
    tmpbuf = IOBuffer()
    blocksize = 0
    higher = ZoomData[]
    compressed = Vector{UInt8}(div(buffer.max_buffer_size * 11, 10))
    dataperblock = div(buffer.max_buffer_size, ZOOM_DATA_SIZE)
    seekstart(buffer.stream1)
    # stream2 is now empty
    @assert position(buffer.stream2) == 0 && eof(buffer.stream2)

    # TODO: Is this correct?
    dataoffset = position(output)
    write(output, UInt32(buffer.count))

    buffer.count = 0
    while !eof(buffer.stream1)
        data = read(buffer.stream1, ZoomData)

        # write data
        write(tmpbuf, data)
        lo = min(lo, (data.chromid, data.chromstart))
        up = max(up, (data.chromid, data.chromend))
        blocksize += 1
        if blocksize == dataperblock
            datasize = compress!(compressed, takebuf_array(tmpbuf))
            push!(blocks, Block(lo, up, position(output), datasize))
            unsafe_write(output, pointer(compressed), datasize)
            lo = (typemax(UInt32), typemax(UInt32))
            up = (typemin(UInt32), typemin(UInt32))
            truncate(tmpbuf, 0)
            blocksize = 0
        end

        # create higher-level zoom
        if length(higher) == scale || (!isempty(higher) && data.chromid != higher[1].chromid)
            write(buffer.stream2, aggregate(higher))
            buffer.count += 1
            empty!(higher)
        end
        push!(higher, data)
    end

    # write remaining data and zoom if any
    if blocksize > 0
        datasize = compress!(compressed, takebuf_array(tmpbuf))
        push!(blocks, Block(lo, up, position(output), datasize))
        unsafe_write(output, pointer(compressed), datasize)
    end
    if !isempty(higher)
        write(buffer.stream2, aggregate(higher))
    end

    # write zoom index
    indexoffset = position(output)
    write_rtree(output, blocks)

    # swap the temporary files
    buffer.stream1, buffer.stream2 = buffer.stream2, buffer.stream1
    seekstart(buffer.stream2)
    truncate(buffer.stream2, 0)

    return dataoffset, indexoffset
end

function aggregate(data::Vector{ZoomData})
    @assert !isempty(data)
    # assumes data are sorted, adjacent, and non-overlapping
    cov = data[1].cov
    min = data[1].min
    max = data[1].max
    sum = data[1].sum
    ssq = data[1].ssq
    for i in 2:endof(data)
        cov += data[i].cov
        min = Base.min(min, data[i].min)
        max = Base.max(max, data[i].max)
        sum += data[i].sum
        ssq += data[i].ssq
    end
    return ZoomData(
        data[1].chromid,
        data[1].chromstart,
        data[end].chromend,
        cov, min, max, sum, ssq)
end
