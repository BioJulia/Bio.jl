# BBI R-Tree
# ==========

const RTREE_MAGIC = 0x2468ACE0

# Supplemental Table 14.
immutable RTreeHeader
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

function Base.write(io::IO, header::RTreeHeader)
    return write(
        io,
        header.magic,
        header.block_size,
        header.item_count,
        header.start_chrom_ix,
        header.start_base,
        header.end_chrom_ix,
        header.end_base,
        header.end_file_offset,
        header.items_per_slot,
        header.reserved)
end

# Supplemental Table 15.
immutable RTreeNode
    isleaf::UInt8
    reserved::UInt8
    count::UInt16
end

function Base.write(io::IO, node::RTreeNode)
    return write(io, node.isleaf, node.reserved, node.count)
end

function isleaf(node::RTreeNode)
    return node.isleaf != 0x00
end

function Base.read(io::IO, ::Type{RTreeNode})
    return RTreeNode(read(io, UInt8), read(io, UInt8), read(io, UInt16))
end

# Supplemental Table 16.
immutable RTreeLeafNode
    start_chrom_ix::UInt32
    start_base::UInt32
    end_chrom_ix::UInt32
    end_base::UInt32
    data_offset::UInt64
    data_size::UInt64
end

# Supplemental Table 17
immutable RTreeInternalNode
    start_chrom_ix::UInt32
    start_base::UInt32
    end_chrom_ix::UInt32
    end_base::UInt32
    data_offset::UInt64
end

# disk-serialized R tree
immutable RTree{T<:IO}
    stream::T
    offset::UInt64
    header::RTreeHeader
end

function RTree(stream::IO, offset::Integer)
    read32() = read(stream, UInt32)
    read64() = read(stream, UInt64)
    seek(stream, offset)
    magic = read32()
    if magic != RTREE_MAGIC
        error("invalid R tree magic bytes")
    end
    return RTree(
        stream,
        offset,
        RTreeHeader(
            magic,    read32(), read64(), read32(), read32(),
            read32(), read32(), read64(), read32(), read32()))
end

function find_overlapping_data_offsets(tree::RTree, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    root_offset = tree.offset + sizeof(RTreeHeader) % UInt64
    stack = [root_offset]
    ret = UInt64[]
    while !isempty(stack)
        offset = pop!(stack)
        seek(tree.stream, offset)
        node = read(tree.stream, RTreeNode)
        for i in 1:node.count
            start_chrom_id = read(tree.stream, UInt32)
            start_chrom_start = read(tree.stream, UInt32)
            end_chrom_id = read(tree.stream, UInt32)
            end_chrom_end = read(tree.stream, UInt32)
            if overlaps(chromid, chromstart, chromend, start_chrom_id, start_chrom_start, end_chrom_id, end_chrom_end)
                offset = read(tree.stream, UInt64)
                if isleaf(node)
                    push!(ret, offset)
                    skip(tree.stream, 8)  # skip data_size
                else
                    push!(stack, offset)
                end
            else
                skip(tree.stream, isleaf(node) ? 16 : 8)
            end
        end
    end
    sort!(ret)
    return ret
end

function overlaps(chromid, chromstart, chromend, start_chrom_id, start_chrom_start, end_chrom_id, end_chrom_end)
    if chromid < start_chrom_id || (chromid == start_chrom_id && chromend ≤ start_chrom_start)
        # query is strictly left
        return false
    elseif chromid > end_chrom_id || (chromid == end_chrom_id && chromstart ≥ end_chrom_end)
        # query is strictly right
        return false
    else
        return true
    end
end


function write_rtree(stream::IO, summaries::Vector{SectionSummary})
    blocksize = 16
    offset = position(stream)
    if isempty(summaries)
        header = RTreeHeader(
            RTREE_MAGIC,
            blocksize,
            0, 0, 0, 0, 0,
            offset, 0, 0)
    else
        header = RTreeHeader(
            RTREE_MAGIC,
            blocksize,
            count_chroms(summaries),
            summaries[1].chromid,
            summaries[1].chromstart,
            summaries[end].chromid,
            summaries[end].chromend,
            offset,
            0,  # how do we set items_per_slot?
            0   # reserved
        )
    end
    # TODO: use internal nodes to reduce data read
    n = 0
    n += write(stream, header)
    node = RTreeNode(0x01, 0x00, length(summaries))
    n += write(stream, node)
    for s in summaries
        n += write(
            stream,
            s.chromid,
            s.chromstart,
            s.chromid,
            s.chromend,
            s.offset,
            # datasize?
            0)
    end
    return n
end

function count_chroms(summaries::Vector{SectionSummary})
    return length(Set(s.chromid for s in summaries))
end
