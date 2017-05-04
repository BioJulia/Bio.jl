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
immutable RTreeNodeFormat
    isleaf::UInt8
    reserved::UInt8
    count::UInt16
end

function Base.write(io::IO, node::RTreeNodeFormat)
    return write(io, node.isleaf, node.reserved, node.count)
end

function isleaf(node::RTreeNodeFormat)
    return node.isleaf != 0x00
end

function Base.read(io::IO, ::Type{RTreeNodeFormat})
    return RTreeNodeFormat(read(io, UInt8), read(io, UInt8), read(io, UInt16))
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

function Base.write(stream::IO, node::RTreeLeafNode)
    return write(
        stream,
        node.start_chrom_ix,
        node.start_base,
        node.end_chrom_ix,
        node.end_base,
        node.data_offset,
        node.data_size)
end

# Supplemental Table 17
immutable RTreeInternalNode
    start_chrom_ix::UInt32
    start_base::UInt32
    end_chrom_ix::UInt32
    end_base::UInt32
    data_offset::UInt64
end

function Base.write(stream::IO, node::RTreeInternalNode)
    return write(
        stream,
        node.start_chrom_ix,
        node.start_base,
        node.end_chrom_ix,
        node.end_base,
        node.data_offset)
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

# Find overlapping leaf nodes.
function find_overlapping_nodes(tree::RTree, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    root_offset = tree.offset + (sizeof(RTreeHeader) % UInt64)
    stack = [root_offset]
    ret = RTreeLeafNode[]
    while !isempty(stack)
        offset = pop!(stack)
        seek(tree.stream, offset)
        node = read(tree.stream, RTreeNodeFormat)
        for i in 1:node.count
            start_chrom_id = read(tree.stream, UInt32)
            start_chrom_start = read(tree.stream, UInt32)
            end_chrom_id = read(tree.stream, UInt32)
            end_chrom_end = read(tree.stream, UInt32)
            if overlaps(chromid, chromstart, chromend, start_chrom_id, start_chrom_start, end_chrom_id, end_chrom_end)
                offset = read(tree.stream, UInt64)
                if isleaf(node)
                    datasize = read(tree.stream, UInt64)
                    push!(
                        ret,
                        RTreeLeafNode(
                            start_chrom_id,
                            start_chrom_start,
                            end_chrom_id,
                            end_chrom_end,
                            offset,
                            datasize))
                else
                    push!(stack, offset)
                end
            else
                skip(tree.stream, isleaf(node) ? 16 : 8)
            end
        end
    end
    # Leaves may be already sorted?
    sort!(ret, by=n->n.data_offset)
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

immutable InMemoryRTree
    start_chromid::UInt32
    start_chromstart::UInt32
    end_chromid::UInt32
    end_chromend::UInt32
    children::Vector{InMemoryRTree}
    offset::UInt64
    datasize::UInt64
end

function isleaf(rtree::InMemoryRTree)
    return rtree.offset == 1
end

function InMemoryRTree(children::Vector{InMemoryRTree})
    @assert !isempty(children)
    chromend = children[1].end_chromend
    for i in 2:endof(children)
        chromend = max(chromend, children[i].end_chromend)
    end
    return InMemoryRTree(
        children[1].start_chromid,
        children[1].start_chromstart,
        children[end].end_chromid,
        chromend,
        children,
        0, 0)
end

function write_rtree(stream::IO, summaries::Vector{SectionSummary})
    blocksize = 64
    items_per_slot = 1
    root = build_inmemory_rtree(summaries, blocksize)
    header = RTreeHeader(
        RTREE_MAGIC,
        blocksize,
        root.end_chromid - root.start_chromid + 1,
        root.start_chromid,
        root.start_chromstart,
        root.end_chromid,
        root.end_chromend,
        position(stream),
        items_per_slot,
        # reserved
        0)
    n = write(stream, header)

    # write R-tree recursively
    function rec(node, offset)
        # compute the size of node
        nodesize = sizeof(RTreeNodeFormat) + (isleaf(node) ? sizeof(RTreeLeafNode) : sizeof(RTreeInternalNode)) * length(node.children)

        # write children
        offsets = UInt64[]
        datasizes = UInt64[]
        child_offset = offset + nodesize
        for child in node.children
            if isleaf(node)
                push!(offsets, child.offset)
                push!(datasizes, child.datasize)
            else
                push!(offsets, child_offset)
                push!(datasizes, UInt64(0))
                child_offset = rec(child, child_offset)
            end
        end

        # write node itself
        truncate(stream, child_offset)  # This is necessary to assert `position(stream) == offset` after seeking.
        seek(stream, offset)
        @assert position(stream) == offset
        write(stream, RTreeNodeFormat(isleaf(node) ? 0x01 : 0x00, 0x00, length(node.children)))
        for i in 1:endof(node.children)
            child = node.children[i]
            if isleaf(node)
                write(
                    stream,
                    RTreeLeafNode(
                        child.start_chromid,
                        child.start_chromstart,
                        child.end_chromid,
                        child.end_chromend,
                        offsets[i],
                        datasizes[i]))
            else
                write(
                    stream,
                    RTreeInternalNode(
                        child.start_chromid,
                        child.start_chromstart,
                        child.end_chromid,
                        child.end_chromend,
                        offsets[i]))
            end
        end

        return child_offset
    end

    root_offset = position(stream)
    return n + (rec(root, root_offset) - root_offset)
end

# Build an in-memory B-Tree from leaves to root (bottom-up).
function build_inmemory_rtree(summaries::Vector{SectionSummary}, blocksize::Int)
    function rec(summaries)
        if length(summaries) ≤ blocksize
            # store indexes in leaves
            children = InMemoryRTree[]
            if isempty(summaries)
                start_chromid = end_chromid = UInt32(0)
                start_chromstart = end_chromend = UInt32(0)
            else
                start_chromid = summaries[1].chromid
                start_chromstart = summaries[1].chromstart
                end_chromid = summaries[end].chromid
                end_chromend = summaries[1].chromend
                for s in summaries
                    end_chromend = max(end_chromend, s.chromend)
                    push!(children, InMemoryRTree(s.chromid, s.chromstart, s.chromid, s.chromend, InMemoryRTree[], s.offset, s.datasize))
                end
            end
            return InMemoryRTree(start_chromid, start_chromstart, end_chromid, end_chromend, children, 1, 0)
        else
            d = cld(length(summaries), blocksize)
            children = InMemoryRTree[]
            for i in 1:blocksize
                idx = (i-1)*d+1:min(i*d,endof(summaries))
                push!(children, rec(view(summaries, idx)))
            end
            return InMemoryRTree(children)
        end
    end
    return rec(view(summaries, 1:endof(summaries)))
end
