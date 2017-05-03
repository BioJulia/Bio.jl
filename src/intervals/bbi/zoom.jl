# BBI Zoom Header/Data
# ====================

# Supplemental Table 6.
immutable ZoomHeader
    reduction_level::UInt32
    reserved::UInt32
    data_offset::UInt64
    index_offset::UInt64
end

function Base.read(io::IO, ::Type{ZoomHeader})
    u32() = read(io, UInt32)
    u64() = read(io, UInt64)
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

function Base.read(io::IO, ::Type{ZoomData})
    u32() = read(io, UInt32)
    f32() = read(io, Float32)
    return ZoomData(
        u32(), u32(), u32(), u32(),
        f32(), f32(), f32(), f32())
end

function find_zoom_data(zoom::Zoom, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    nodes = find_overlapping_nodes(zoom.rtree, chromid, chromstart, chromend)
    stream = zoom.rtree.stream
    ret = ZoomData[]
    for node in nodes
        seek(stream, node.data_offset)
        while position(stream) < node.data_offset + node.data_size
            data = read(stream, ZoomData)
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

function compare_intervals(data::ZoomData, interval::Tuple{UInt32,UInt32,UInt32})
    # TODO
end
