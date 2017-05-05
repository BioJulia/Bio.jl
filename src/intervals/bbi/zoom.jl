# BBI Zoom Header/Data
# ====================

# Supplemental Table 6.
immutable ZoomHeader
    reduction_level::UInt32
    reserved::UInt32
    data_offset::UInt64
    index_offset::UInt64
end

function Base.read(stream::IO, ::Type{ZoomHeader})
    u32() = read(stream, UInt32)
    u64() = read(stream, UInt64)
    return ZoomHeader(u32(), u32(), u64(), u64())
end

immutable Zoom{T<:IO}
    header::ZoomHeader
    rtree::RTree{T}
end

function Zoom(stream::IO, header::ZoomHeader)
    return Zoom(header, RTree(stream, header.index_offset))
end

# Supplemental Table 19.
immutable ZoomData
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
    count::UInt32
    minval::Float32
    maxval::Float32
    sumval::Float32
    sumsqval::Float32
end

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
        data.chromid,
        data.chromstart,
        data.chromend,
        data.count,
        data.minval,
        data.maxval,
        data.sumval,
        data.sumsqval)
end

function find_overlapping_zoomdata(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    nodes = find_overlapping_nodes(zoom.rtree, chromid, chromstart, chromend)
    stream = zoom.rtree.stream
    ret = ZoomData[]
    for node in nodes
        seek(stream, node.offset)
        # TODO: lazy decompression
        datastream = IOBuffer(Libz.decompress(read(stream, node.size)))
        while !eof(datastream)
            data = read(datastream, ZoomData)
            c = compare_intervals(data, (chromid, chromstart, chromend))
            if c == 0  # overlap
                push!(ret, data)
            elseif c > 0
                break
            end
        end
    end
    return ret
end

function compare_intervals(data::ZoomData, interval::NTuple{3,UInt32})
    chromid, chromstart, chromend = interval
    if data.chromid < chromid || (data.chromid == chromid && data.chromend ≤ chromstart)
        # strictly left
        return -1
    elseif data.chromid > chromid || (data.chromid == chromid && data.chromstart ≥ chromend)
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
        count += round(UInt32, d.count * cov / (d.chromend - d.chromstart))
    end
    return count
end

function mean(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    sum = 0.0f0
    size = 0
    for d in data
        cov = coverage2((d.chromstart, d.chromend), (chromstart, chromend))
        sum += (d.sumval * cov) / (d.chromend - d.chromstart)
        size += cov
    end
    return sum / size
end

function minimum(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    return isempty(data) ? NaN32 : foldl((x, d) -> min(x, d.minval), +Inf32, data)
end

function maximum(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    return isempty(data) ? NaN32 : foldl((x, d) -> max(x, d.maxval), -Inf32, data)
end

function std(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    sum = 0.0f0
    ssq = 0.0f0
    size = 0
    for d in data
        cov = coverage2((d.chromstart, d.chromend), (chromstart, chromend))
        sum += (d.sumval * cov) / (d.chromend - d.chromstart)
        ssq += (d.sumsqval * cov) / (d.chromend - d.chromstart)
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


# Zoom Buffer
# -----------

type ZoomBuffer
    # output stream
    stream::IO

    # chrom ID => length
    chromlens::Dict{UInt32,UInt32}

    # bin size of zoom
    binsize::UInt32

    # number of data per block
    dataperblock::Int

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

    # flag whether the zoom statistics are written to the output stream
    written::Bool
end

function ZoomBuffer(stream, chromlist, binsize, uncompressed_buffer_size)
    return ZoomBuffer(
        stream,
        Dict(id => len for (name, id, len) in chromlist),
        binsize,
        fld(uncompressed_buffer_size, sizeof(ZoomData)),
        typemax(UInt32),
        Block[],
        UInt32[], Float32[], Float32[], Float32[], Float32[],
        true)
end

function add_value!(buffer::ZoomBuffer, chromid::UInt32, chromstart::UInt32, chromend::UInt32, value::Float32)
    if buffer.chromid == typemax(UInt32)
        init_buffer!(buffer, chromid)
    elseif chromid != buffer.chromid
        write_buffer(buffer)
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
    buffer.written = false
    return buffer
end

function init_buffer!(buffer::ZoomBuffer, chromid::UInt32)
    @assert buffer.written == true
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

function write_buffer(buffer::ZoomBuffer)
    @assert buffer.chromid != typemax(UInt32)
    # write zoom data
    nbins = length(buffer.cov)
    binsize = buffer.binsize
    dataperblock = buffer.dataperblock
    for i in 1:cld(nbins, dataperblock)
        blockbuf = IOBuffer()
        binstart = (i-1)*dataperblock
        binend = min(i*dataperblock, nbins)
        for bin in binstart:binend-1
            chromstart = bin * binsize
            bin1 = bin + 1
            write(
                blockbuf,
                ZoomData(
                    buffer.chromid,
                    chromstart,
                    chromstart + binsize,
                    buffer.cov[bin1],
                    buffer.min[bin1],
                    buffer.max[bin1],
                    buffer.sum[bin1],
                    buffer.ssq[bin1]))
        end
        lo = (buffer.chromid, binstart * binsize)
        up = (buffer.chromid, binend * binsize)
        data = Libz.compress(takebuf_array(blockbuf))
        push!(buffer.blocks, Block(lo, up, position(buffer.stream), sizeof(data)))
        write(buffer.stream, data)
    end
    buffer.written = true
    return
end

# Write zoom data to `output` from `buffer`.
function write_zoom(output::IO, buffer::ZoomBuffer)
    if !buffer.written
        write_buffer(buffer)
    end
    @assert buffer.written

    # adjust offsets in blocks to output
    offset = position(output)
    blocks′ = similar(buffer.blocks)
    for i in 1:endof(buffer.blocks)
        block = buffer.blocks[i]
        blocks′[i] = Block(block.lo, block.up, block.offset + offset, block.size)
    end

    # write zoom data
    seekstart(buffer.stream)
    write(output, buffer.stream)

    # write R-tree index
    write_rtree(output, blocks′)

    return
end
