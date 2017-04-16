# BBI B+ Tree
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

# Supplemental Table 9.
immutable BTreeNode
    isleaf::UInt8
    reserved::UInt8
    count::UInt16
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
    return BTree(stream, offset, BTreeHeader(magic, read32(), read32(), read32(), read64(), read64()))
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
