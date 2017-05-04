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

# Get regional statistics using the zoom data.
function zoomstats(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    data = find_overlapping_zoomdata(zoom, chromid, chromstart, chromend)
    if isempty(data)
        return ZoomData(chromid, chromstart, chromend, 0, NaN32, NaN32, NaN32, NaN32)
    end
    count = UInt32(0)
    minval = Inf32
    maxval = -Inf32
    sumval = 0.0f0
    sumsqval = 0.0f0
    for d in data
        # TODO: should this be normalized by effective length?
        count += d.count
        minval = min(minval, d.minval)
        maxval = max(maxval, d.maxval)
        sumval += d.sumval
        sumsqval += d.sumsqval
    end
    return ZoomData(chromid, chromstart, chromend, count, minval, maxval, sumval, sumsqval)
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
