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

function find_overlapping_zoomdata(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    nodes = find_overlapping_nodes(zoom.rtree, chromid, chromstart, chromend)
    stream = zoom.rtree.stream
    ret = ZoomData[]
    for node in nodes
        seek(stream, node.data_offset)
        # TODO: lazy decompression
        datastream = IOBuffer(Libz.decompress(read(stream, node.data_size)))
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
