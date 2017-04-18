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
    chrom_id::UInt32
    chrom_start::UInt32
    chrom_end::UInt32
    valid_count::UInt32
    min::Float32
    max::Float32
    sum::Float32
    ssq::Float32
end

function Base.read(io::IO, ::Type{ZoomData})
    u32() = read(io, UInt32)
    f32() = read(io, Float32)
    return ZoomData(
        u32(), u32(), u32(), u32(),
        f32(), f32(), f32(), f32())
end
