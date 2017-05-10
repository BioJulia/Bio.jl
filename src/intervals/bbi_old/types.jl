# Types
# =====

const BIGBED_MAGIC = 0x8789F2EB
const BIGWIG_MAGIC = 0x888FFC26
const BIGBED_BTREE_MAGIC = 0x78CA8C91
const BIGBED_RTREE_MAGIC = 0x2468ACE0

const BIGBED_MAX_ZOOM_LEVELS = 10
const BIGBED_CURRENT_VERSION = 4

const BIGWIG_DATATYPE_BEDGRAPH  = UInt8(1)
const BIGWIG_DATATYPE_VARSTEP   = UInt8(2)
const BIGWIG_DATATYPE_FIXEDSTEP = UInt8(3)

# See Supplemental Table 5
immutable BigBedHeader
    magic::UInt32
    version::UInt16
    zoom_levels::UInt16
    chromosome_tree_offset::UInt64
    full_data_offset::UInt64
    full_index_offset::UInt64
    field_count::UInt16
    defined_field_count::UInt16
    auto_sql_offset::UInt64
    total_summary_offset::UInt64
    uncompress_buf_size::UInt32
    reserved::UInt64
end

function Base.read(io::IO, ::Type{BigBedHeader})
    return BigBedHeader(
        read(io, UInt32), read(io, UInt16), read(io, UInt16),
        read(io, UInt64), read(io, UInt64), read(io, UInt64),
        read(io, UInt16), read(io, UInt16), read(io, UInt64),
        read(io, UInt64), read(io, UInt32), read(io, UInt64))
end

function Base.write(io::IO, header::BigBedHeader)
    write(io, header.magic)
    write(io, header.version)
    write(io, header.zoom_levels)
    write(io, header.chromosome_tree_offset)
    write(io, header.full_data_offset)
    write(io, header.full_index_offset)
    write(io, header.field_count)
    write(io, header.defined_field_count)
    write(io, header.auto_sql_offset)
    write(io, header.total_summary_offset)
    write(io, header.uncompress_buf_size)
    write(io, header.reserved)
end

# See Supplemental Table 6
immutable BigBedZoomHeader
    reduction_level::UInt32
    reserved::UInt32
    data_offset::UInt64
    index_offset::UInt64
end

function Base.read(io::IO, ::Type{BigBedZoomHeader})
    return BigBedZoomHeader(
        read(io, UInt32), read(io, UInt32), read(io, UInt64), read(io, UInt64))
end

function Base.write(io::IO, header::BigBedZoomHeader)
    write(io, header.reduction_level)
    write(io, header.reserved)
    write(io, header.data_offset)
    write(io, header.index_offset)
end

# Supplemental Table 7
immutable BigBedTotalSummary
    bases_covered::UInt64
    min_val::Float64
    max_val::Float64
    sum_data::Float64
    sum_squares::Float64
end

function Base.read(io::IO, ::Type{BigBedTotalSummary})
    return BigBedTotalSummary(
        read(io, UInt64), read(io, UInt64), read(io, UInt64),
        read(io, UInt64), read(io, UInt64))
end

function Base.write(io::IO, summary::BigBedTotalSummary)
    write(io, summary.bases_covered)
    write(io, summary.min_val)
    write(io, summary.max_val)
    write(io, summary.sum_data)
    write(io, summary.sum_squares)
end

# Supplemental Table 19
immutable BigBedZoomData
    chrom_id::UInt32
    chrom_start::UInt32
    chrom_end::UInt32
    valid_count::UInt32
    min_val::Float32
    max_val::Float32
    sum_data::Float32
    sum_squares::Float32
end

function Base.read(io::IO, ::Type{BigBedZoomData})
    return BigBedZoomData(
        read(io, UInt32), read(io, UInt32), read(io, UInt32), read(io, UInt32),
        read(io, Float32), read(io, Float32), read(Float32), read(Float32))
end

function Base.write(io::IO, data::BigBedZoomData)
    write(io, data.chrom_id)
    write(io, data.chrom_start)
    write(io, data.chrom_end)
    write(io, data.valid_count)
    write(io, data.min_val)
    write(io, data.max_val)
    write(io, data.sum_data)
    write(io, data.sum_squares)
end

# Supplemental Table 8
immutable BigBedBTreeHeader
    magic::UInt32
    block_size::UInt32
    key_size::UInt32
    val_size::UInt32
    item_count::UInt64
    reserved::UInt64
end

function Base.read(io::IO, ::Type{BigBedBTreeHeader})
    return BigBedBTreeHeader(
        read(io, UInt32), read(io, UInt32), read(io, UInt32),
        read(io, UInt32), read(io, UInt64), read(io, UInt64))
end

function Base.write(io::IO, header::BigBedBTreeHeader)
    write(io, header.magic, header.block_size, header.key_size,
          header.val_size, header.item_count, header.reserved)
end

# Supplementary Table 9
immutable BigBedBTreeNode
    isleaf::UInt8
    reserved::UInt8
    count::UInt16
end

function Base.read(io::IO, ::Type{BigBedBTreeNode})
    return BigBedBTreeNode(read(io, UInt8), read(io, UInt8), read(io, UInt16))
end

function Base.write(io::IO, node::BigBedBTreeNode)
    write(io, node.isleaf)
    write(io, node.reserved)
    write(io, node.count)
end

# Supplemental Table 10
type BigBedBTreeLeafNode
    key::Vector{UInt8}
    chrom_id::UInt32
    chrom_size::UInt32
end

function BigBedBTreeLeafNode(keysize::Integer)
    BigBedBTreeLeafNode(Array(UInt8, keysize), 0, 0)
end

function Base.read!(io::IO, node::BigBedBTreeLeafNode)
    nb = readbytes!(io, node.key, length(node.key))
    if nb < length(node.key)
        error("Unexpected end of input.")
    end
    node.chrom_id = read(io, UInt32)
    node.chrom_size = read(io, UInt32)
    return nb + 8
end

# Supplemental Table 11
type BigBedBTreeInternalNode
    key::Vector{UInt8}
    child_offset::UInt64
end

function BigBedBTreeInternalNode(keysize::Integer)
    BigBedBTreeInternalNode(Array(UInt8, keysize), 0)
end

function Base.read!(io::IO, node::BigBedBTreeInternalNode)
    if readbytes!(io, node.key, length(node.key)) < length(node.key)
        error("Unexpected end of input.")
    end
    node.child_offset = read(io, UInt64)
end

# Supplemental Table 13
immutable BigWigSectionHeader
    chrom_id::UInt32
    chrom_start::UInt32
    chrom_end::UInt32
    item_step::UInt32
    item_span::UInt32
    data_type::UInt8
    reserved::UInt8
    item_count::UInt16
end

function Base.write(io::IO, header::BigWigSectionHeader)
    write(io, header.chrom_id)
    write(io, header.chrom_start)
    write(io, header.chrom_end)
    write(io, header.item_step)
    write(io, header.item_span)
    write(io, header.data_type)
    write(io, header.reserved)
    write(io, header.item_count)
end

immutable BedGraphItem
    chrom_start::UInt32
    chrom_end::UInt32
    val::Float32
end

function Base.write(io::IO, item::BedGraphItem)
    write(io, item.chrom_start)
    write(io, item.chrom_end)
    write(io, item.val)
end

# Supplemental Table 14
immutable BigBedRTreeHeader
    magic::UInt32
    block_size::UInt32
    item_count::UInt64
    start_chrom_ix::UInt32
    start_base::UInt32
    end_chrom_ix::UInt32
    end_base::UInt32
    end_file_offset::UInt64
    items_per_slot::UInt32
    reserved::UInt32
end

function Base.read(io::IO, ::Type{BigBedRTreeHeader})
    return BigBedRTreeHeader(
        read(io, UInt32), read(io, UInt32), read(io, UInt64), read(io, UInt32),
        read(io, UInt32), read(io, UInt32), read(io, UInt32), read(io, UInt64),
        read(io, UInt32), read(io, UInt32))
end

function Base.write(io::IO, header::BigBedRTreeHeader)
    write(io, header.magic)
    write(io, header.block_size)
    write(io, header.item_count)
    write(io, header.start_chrom_ix)
    write(io, header.start_base)
    write(io, header.end_chrom_ix)
    write(io, header.end_base)
    write(io, header.end_file_offset)
    write(io, header.items_per_slot)
    write(io, header.reserved)
end

# Supplemental Table 15
immutable BigBedRTreeNode
    isleaf::UInt8
    reserved::UInt8
    count::UInt16
end

function Base.read(io::IO, ::Type{BigBedRTreeNode})
    return BigBedRTreeNode(read(io, UInt8), read(io, UInt8), read(io, UInt16))
end

function Base.write(io::IO, node::BigBedRTreeNode)
    write(io, node.isleaf)
    write(io, node.reserved)
    write(io, node.count)
end

# Supplemental Table 16
immutable BigBedRTreeLeafNode
    start_chrom_ix::UInt32
    start_base::UInt32
    end_chrom_ix::UInt32
    end_base::UInt32
    data_offset::UInt64
    data_size::UInt64
end

function Base.read(io::IO, ::Type{BigBedRTreeLeafNode})
    return BigBedRTreeLeafNode(
        read(io, UInt32), read(io, UInt32), read(io, UInt32), read(io, UInt32),
        read(io, UInt64), read(io, UInt64))
end

function Base.write(io::IO, node::BigBedRTreeLeafNode)
    write(io, node.start_chrom_ix)
    write(io, node.start_base)
    write(io, node.end_chrom_ix)
    write(io, node.end_base)
    write(io, node.data_offset)
    write(io, node.data_size)
end

# Supplemental Table 17
immutable BigBedRTreeInternalNode
    start_chrom_ix::UInt32
    start_base::UInt32
    end_chrom_ix::UInt32
    end_base::UInt32
    data_offset::UInt64
end

function Base.read(io::IO, ::Type{BigBedRTreeInternalNode})
    return BigBedRTreeInternalNode(
        read(io, UInt32), read(io, UInt32), read(io, UInt32), read(io, UInt32),
        read(io, UInt64))
end

function Base.write(io::IO, node::BigBedRTreeInternalNode)
    write(io, node.start_chrom_ix)
    write(io, node.start_base)
    write(io, node.end_chrom_ix)
    write(io, node.end_base)
    write(io, node.data_offset)
end
