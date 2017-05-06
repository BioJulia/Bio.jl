# BBI R-Tree
# ==========

const RTREE_MAGIC = 0x2468ACE0

# Supplemental Table 14.
immutable RTreeHeader
    magic::UInt32
    block_size::UInt32
    item_count::UInt64
    lo::NTuple{2,UInt32}
    up::NTuple{2,UInt32}
    end_file_offset::UInt64
    items_per_slot::UInt32
    reserved::UInt32
end

const RTREE_HEADER_SIZE = 48

function Base.read(stream::IO, ::Type{RTreeHeader})
    return RTreeHeader(
        read(stream, UInt32),
        read(stream, UInt32),
        read(stream, UInt64),
        readbound(stream),
        readbound(stream),
        read(stream, UInt64),
        read(stream, UInt32),
        read(stream, UInt32))
end

function Base.write(stream::IO, header::RTreeHeader)
    return write(
        stream,
        header.magic,
        header.block_size,
        header.item_count,
        header.lo[1], header.lo[2],
        header.up[1], header.up[2],
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

const RTREE_NODE_FORMAT_SIZE = 4

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
    lo::NTuple{2,UInt32}
    up::NTuple{2,UInt32}
    offset::UInt64
    size::UInt64
end

const RTREE_LEAF_NODE_SIZE = 32

function Base.write(stream::IO, node::RTreeLeafNode)
    return write(
        stream,
        node.lo[1], node.lo[2],
        node.up[1], node.up[2],
        node.offset, node.size)
end

# Supplemental Table 17
immutable RTreeInternalNode
    lo::NTuple{2,UInt32}
    up::NTuple{2,UInt32}
    offset::UInt64
end

const RTREE_INTERNAL_NODE_SIZE = 32

function Base.write(stream::IO, node::RTreeInternalNode)
    return write(
        stream,
        node.lo[1], node.lo[2],
        node.up[1], node.up[2],
        node.offset)
end


# TODO: Merge this with RTreeLeafNode and SectionSummary
# Data block indexed by R-tree (0-origin, left-closed and right-open).
immutable Block
    # lower bound of genomic interval (chromid, chromstart)
    lo::Tuple{UInt32,UInt32}

    # upper bound of genomic interval (chromid, chromend)
    up::Tuple{UInt32,UInt32}

    # offset to data block
    offset::UInt64

    # size of data block (likely compressed, disk dize)
    size::UInt64
end

# disk-serialized R tree
immutable RTree{T<:IO}
    stream::T
    offset::UInt64
    header::RTreeHeader
end

function RTree(stream::IO, offset::Integer)
    seek(stream, offset)
    header = read(stream, RTreeHeader)
    if header.magic != RTREE_MAGIC
        error("invalid R tree magic bytes")
    end
    return RTree(stream, offset, header)
end

# Find overlapping blocks.
function find_overlapping_blocks(tree::RTree, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    root_offset = tree.offset + (RTREE_HEADER_SIZE % UInt64)
    stack = [root_offset]
    ret = Block[]
    while !isempty(stack)
        offset = pop!(stack)
        seek(tree.stream, offset)
        node = read(tree.stream, RTreeNodeFormat)
        for i in 1:node.count
            lo = readbound(tree.stream)
            up = readbound(tree.stream)
            if overlaps(chromid, chromstart, chromend, lo, up)
                offset = read(tree.stream, UInt64)
                if isleaf(node)
                    datasize = read(tree.stream, UInt64)
                    push!(ret, Block(lo, up, offset, datasize))
                else
                    push!(stack, offset)
                end
            else
                skip(tree.stream, isleaf(node) ? 16 : 8)
            end
        end
    end
    # Leaves may be already sorted?
    sort!(ret, by=n->n.offset)
    return ret
end

function overlaps(chromid, chromstart, chromend, lo, up)
    if chromid < lo[1] || (chromid == lo[1] && chromend ≤ lo[2])
        # query is strictly left
        return false
    elseif chromid > up[1] || (chromid == up[1] && chromstart ≥ up[2])
        # query is strictly right
        return false
    else
        return true
    end
end

immutable InMemoryRTree
    lo::NTuple{2,UInt32}
    up::NTuple{2,UInt32}
    offset::UInt64
    size::UInt64
    children::Vector{InMemoryRTree}

    function InMemoryRTree(lo, up, offset, size, children=nothing)
        if children == nothing
            return new(lo, up, offset, size)
        else
            return new(lo, up, offset, size, children)
        end
    end
end

function isleaf(rtree::InMemoryRTree)
    return rtree.offset == 1
end

function InMemoryRTree(children::Vector{InMemoryRTree})
    @assert !isempty(children)
    # TODO: assert this
    #@assert issorted(n.lo for n in children)
    return InMemoryRTree(
        #children[1].lo,
        Base.minimum(n.lo for n in children),
        Base.maximum(n.up for n in children),
        0, 0,
        children)
end

function write_rtree(stream::IO, blocks::Vector{Block})
    blocksize = 64
    items_per_slot = 1
    root = build_inmemory_rtree(blocks, blocksize)
    header = RTreeHeader(
        RTREE_MAGIC,
        blocksize,
        root.up[1] - root.lo[1] + 1,
        root.lo,
        root.up,
        position(stream),
        items_per_slot,
        # reserved
        0)
    n = write(stream, header)

    # write R-tree recursively
    function rec(node, offset)
        # compute the size of node
        nodesize = RTREE_NODE_FORMAT_SIZE + (isleaf(node) ? RTREE_LEAF_NODE_SIZE : RTREE_INTERNAL_NODE_SIZE) * length(node.children)

        # write children
        offsets = UInt64[]
        datasizes = UInt64[]
        child_offset = offset + nodesize
        for child in node.children
            if isleaf(node)
                push!(offsets, child.offset)
                push!(datasizes, child.size)
            else
                push!(offsets, child_offset)
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
                write(stream, RTreeLeafNode(child.lo, child.up, offsets[i], datasizes[i]))
            else
                write(stream, RTreeInternalNode(child.lo, child.up, offsets[i]))
            end
        end

        return child_offset
    end

    root_offset = position(stream)
    return n + (rec(root, root_offset) - root_offset)
end

# Build an in-memory B-Tree from leaves to root (bottom-up).
function build_inmemory_rtree(blocks::Vector{Block}, blocksize::Int)
    if !issorted(blocks, by=b->b.lo)
        blocks = sort(blocks, by=b->b.lo)
    end
    function rec(blocks)
        if length(blocks) ≤ blocksize
            # store indexes in leaves
            children = InMemoryRTree[]
            if isempty(blocks)
                lo = up = (UInt32(0), UInt32(0))
            else
                lo = blocks[1].lo
                up = blocks[1].up
                for block in blocks
                    up = max(up, block.up)
                    push!(children, InMemoryRTree(block.lo, block.up, block.offset, block.size))
                end
            end
            return InMemoryRTree(lo, up, 1, 0, children)
        else
            d = cld(length(blocks), blocksize)
            children = InMemoryRTree[]
            for i in 1:blocksize
                idx = (i-1)*d+1:min(i*d,endof(blocks))
                push!(children, rec(view(blocks, idx)))
            end
            return InMemoryRTree(children)
        end
    end
    return rec(view(blocks, 1:endof(blocks)))
end

function readbound(stream::IO)
    return read(stream, UInt32), read(stream, UInt32)
end
