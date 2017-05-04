# BBI B+-Tree
# ===========

const BTREE_MAGIC = 0x78CA8C91

# Supplemental Table 8.
immutable BTreeHeader
    magic::UInt32
    block_size::UInt32
    key_size::UInt32
    val_size::UInt32
    item_count::UInt64
    reserved::UInt64
end

function Base.write(stream::IO, header::BTreeHeader)
    return write(
        stream,
        header.magic,    header.block_size, header.key_size,
        header.val_size, header.item_count, header.reserved)
end

# Supplemental Table 9.
immutable BTreeNode
    isleaf::UInt8
    reserved::UInt8
    count::UInt16
end

function Base.write(stream::IO, node::BTreeNode)
    return write(stream, node.isleaf, node.reserved, node.count)
end

function isleaf(node::BTreeNode)
    return node.isleaf != 0
end

function Base.read(io::IO, ::Type{BTreeNode})
    return BTreeNode(read(io, UInt8), read(io, UInt8), read(io, UInt16))
end

# disk-serialized B+ tree
immutable BTree{T<:IO}
    stream::T
    offset::UInt64
    header::BTreeHeader
end

function BTree(stream::IO, offset::Integer)
    read32() = read(stream, UInt32)
    read64() = read(stream, UInt64)
    seek(stream, offset)
    magic = read32()
    if magic != BTREE_MAGIC
        error("invalid B+-tree magic bytes")
    end
    return BTree(stream, UInt64(offset), BTreeHeader(magic, read32(), read32(), read32(), read64(), read64()))
end

# Load the chromosome list.
function chromlist(tree::BTree)
    list = Tuple{String,UInt32,UInt32}[]
    # traverse tree
    key = Vector{UInt8}(tree.header.key_size)
    stack = UInt64[tree.offset + sizeof(BTreeHeader)]
    while !isempty(stack)
        offset = pop!(stack)
        seek(tree.stream, offset)
        node = read(tree.stream, BTreeNode)
        for i in 1:node.count
            if isleaf(node)
                read!(tree.stream, key)
                i = findfirst(key, 0x00)
                chromname = i == 0 ? String(copy(key)) : String(key[1:i-1])
                chromid = read(tree.stream, UInt32)
                chromsize = read(tree.stream, UInt32)
                push!(list, (chromname, chromid, chromsize))
            else
                skip(tree.stream, tree.header.key_size)
                offset = read(tree.stream, UInt64)
                push!(stack, offset)
            end
        end
    end
    return list
end

function write_btree(stream::IO, chromlist::Vector{Tuple{String,UInt32,UInt32}})
    # This function stores all chromosomes in the root node as a leaf because it
    # is simple to implement.
    blksize = length(chromlist)
    keysize = Base.maximum(sizeof(name) for (name, _) in chromlist)
    valsize = 8
    n = 0

    # write header
    n += write(stream, BTreeHeader(BTREE_MAGIC, blksize, keysize, valsize, length(chromlist), 0))

    # write the root node format
    n += write(stream, BTreeNode(0x01, 0x00, length(chromlist)))

    # write the root node
    key = Vector{UInt8}(keysize)
    for (name, id, len) in sort(chromlist, by=x->x[1])  # sort by name
        fill!(key, 0x00)
        @assert sizeof(name) ≤ keysize
        Mem.copy(key, name, sizeof(name))
        n += write(stream, key, id, len)
    end

    return n
end

# Add a unique ID for each chromosome; chromlist is a tuple of (chrom name, chrom length).
function add_chrom_ids(chromlist::Vector{Tuple{String,UInt32}})
    chromlist = sort(chromlist, by=first)  # sort by chromosome name
    return [(name, UInt32(id - 1), UInt32(len)) for (id, (name, len)) in enumerate(chromlist)]
end

function add_chrom_ids(chromlist)
    return add_chrom_ids([(String(name), UInt32(len)) for (name, len) in chromlist])
end
