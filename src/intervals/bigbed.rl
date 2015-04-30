

const BIGBED_MAGIC = 0x8789F2EB
const BIGWIG_MAGIC = 0x888FFC26


# See Supplemental Table 5
immutable BigBedHeader
    magic::Uint32
    version::Uint16
    zoom_levels::Uint16
    chromosome_tree_offset::Uint64
    full_data_offset::Uint64
    full_index_offset::Uint64
    field_count::Uint16
    defined_field_count::Uint16
    auto_sql_offset::Uint64
    total_summary_offset::Uint64
    uncompressed_buf_size::Uint32
    reserved::Uint64
end


function read(io::IO, ::Type{BigBedHeader})
    return BigBedHeader(
        read(io, Uint32), read(io, Uint16), read(io, Uint16),
        read(io, Uint64), read(io, Uint64), read(io, Uint64),
        read(io, Uint16), read(io, Uint16), read(io, Uint64),
        read(io, Uint64), read(io, Uint32), read(io, Uint64))
end


# See Supplemental Table 6
immutable BigBedZoomHeader
    reduction_level::Uint32
    reserved::Uint32
    data_offset::Uint64
    index_offset::Uint64
end


function read(io::IO, ::Type{BigBedZoomHeader})
    return BigBedZoomHeader(
        read(io, Uint32), read(io, Uint32), read(io, Uint64), read(io, Uint64))
end


# Supplemental Table 7
immutable BigBedTotalSummary
    bases_covered::Uint64
    minval::Float64
    maxval::Float64
    sumdata::Float64
    sumsquares::Float64
end


function read(io::IO, ::Type{BigBedTotalSummary})
    return BigBedTotalSummary(
        read(io, Uint64), read(io, Uint64), read(io, Uint64),
        read(io, Uint64), read(io, Uint64))
end


# Supplemental Table 8
immutable BigBedBTreeHeader
    magic::Uint32
    block_size::Uint32
    key_size::Uint32
    val_size::Uint32
    item_count::Uint64
    reserved::Uint64
end


function read(io::IO, ::Type{BigBedBTreeHeader})
    return BigBedBTreeHeader(
        read(io, Uint32), read(io, Uint32), read(io, Uint32),
        read(io, Uint32), read(io, Uint64), read(io, Uint64))
end


# Supplementary Table 9
immutable BigBedBTreeNode
    isleaf::Uint8
    reserved::Uint8
    count::Uint16
end


function read(io::IO, ::Type{BigBedBTreeNode})
    return BigBedBTreeNode(read(io, Uint8), read(io, Uint8), read(io, Uint16))
end


# Supplemental Table 10
type BigBedBTreeLeafNode
    key::Vector{Uint8}
    chrom_id::Uint32
    chrom_size::Uint32
end


function BigBedBTreeLeafNode(keysize::Integer)
    BigBedBTreeLeafNode(Array(Uint8, keysize), 0, 0)
end


function read!(io::IO, node::BigBedBTreeLeafNode)
    if readbytes!(io, node.key, length(node.key)) < length(node.key)
        error("Unexpected end of input.")
    end
    node.chrom_id = read(io, Uint32)
    node.chrom_size = read(io, Uint32)
end


# Supplemental Table 11
type BigBedBTreeInternalNode
    key::Vector{Uint8}
    child_offset::Uint64
end


function BigBedBTreeInternalNode(keysize::Integer)
    BigBedBTreeInternalNode(Array(Uint8, keysize), 0)
end


function read!(io::IO, node::BigBedBTreeInternalNode)
    if readbytes!(io, node.key, length(node.key)) < length(node.key)
        error("Unexpected end of input.")
    end
    node.child_offset = read(io, Uint64)
end



using Bio: BufferedReader

immutable BigBed <: FileFormat end


type BigBedData
    reader::BufferedReader
    header::BigBedHeader
    zoom_headers::Vector{BigBedZoomHeader}
    autosql::String
    summary::BigBedTotalSummary
    chromosome_tree_header::BigBedBTreeHeader
    data_count::Uint32

    # preallocated space for reading and searchig the B-tree
    btree_internal_nodes::Vector{BigBedBTreeInternalNode}
    btree_leaf_nodes::Vector{BigBedBTreeLeafNode}
    key::Vector{Uint8}
    node_keys::Vector{Vector{Uint8}}
end


@doc """
Open a BigBed file for reading.

Once opened, entries can be read from the file either by iterating over it, or
by indexing into it with an interval.
""" ->
function read(input::Union(IO, String, Vector{Uint8}),
              ::Type{BigBed}; memory_map::Bool=false)
    reader = BufferedReader(input, memory_map)

    if isa(input, IO) && !applicable(seek, input)
        error("BigBed files can only be read from seekable input, such as regular files.")
    end

    # header
    header = read(reader, BigBedHeader)
    if header.magic != BIGBED_MAGIC
        error("Input is not a valid BigBed file.")
    end

    if header.version < 3
        error("Input is an older unsupported version of BigBed.")
    end

    # zoom headers
    zoom_headers = Array(BigBedZoomHeader, header.zoom_levels)
    for i in 1:header.zoom_levels
        zoom_headers[i] = read(reader, BigBedZoomHeader)
    end

    # autosql
    seek(reader, header.auto_sql_offset + 1)
    autosql_buf = IOBuffer()
    while (c = read(reader, Uint8)) != 0x00
        write(autosql_buf, c)
    end
    # TODO: eventually we should parse this and do something useful with it
    autosql = takebuf_string(autosql_buf)

    # total summary
    seek(reader, header.total_summary_offset + 1)
    summary = read(reader, BigBedTotalSummary)

    # b-tree header
    seek(reader, header.chromosome_tree_offset + 1)
    chromosome_tree_header = read(reader, BigBedBTreeHeader)

    # skip over chromosome tree
    seek(reader, header.full_data_offset + 1)

    data_count = read(reader, Uint32)

    return BigBedData(reader, header, zoom_headers, autosql, summary,
                      chromosome_tree_header, data_count,
                      BigBedBTreeInternalNode[BigBedBTreeInternalNode(chromosome_tree_header.key_size)
                                              for _ in 1:chromosome_tree_header.block_size],
                      BigBedBTreeLeafNode[BigBedBTreeLeafNode(chromosome_tree_header.key_size)
                                          for _ in 1:chromosome_tree_header.block_size],
                      Array(Uint8, chromosome_tree_header.key_size),
                      Array(Vector{Uint8}, chromosome_tree_header.block_size))
end


function memisless(a::Vector{Uint8}, b::Vector{Uint8})
    if length(a) != length(b)
        return length(a) < length(b)
    end
    i = 1
    while i <= length(a)
        if a[i] != b[i]
            return a[i] < b[i]
        end
        i += 1
    end
    return false
end


@doc """
An iterator over all entries in a BigBed file.
""" ->
immutable BigBedIterator
end

# TODO: Linear iterator over all BigBed elements.


@doc """
Find the given seqname in the BigBed file's index and read the corresponding 
sequence id and length
""" ->
function lookup_seqname(bb::BigBedData, seqname::String)
    seek(bb.reader, bb.header.chromosome_tree_offset + sizeof(BigBedBTreeHeader) + 1)

    fill!(bb.key, 0)
    copy!(bb.key, seqname)

    while true
        node = read(bb.reader, BigBedBTreeNode)
        if node.isleaf == 0
            for i in 1:bb.chromosome_tree_header.block_size
                read!(bb.reader, bb.btree_internal_nodes[i])
                bb.node_keys[i] = bb.btree_internal_nodes[i].key
            end
            i = searchsortedfirst(bb.node_keys, bb.key, Base.Order.Lt(memisless))
            if !(1 <= i <= bb.chromosome_tree_header.block_size)
                break
            end

            seek(bb.reader, bb.btree_internal_nodes[i].child_offset)
        else
            for i in 1:bb.chromosome_tree_header.block_size
                read!(bb.reader, bb.btree_leaf_nodes[i])
                bb.node_keys[i] = bb.btree_leaf_nodes[i].key
            end
            i = searchsortedfirst(bb.node_keys, bb.key, Base.Order.Lt(memisless))
            if !(1 <= i <= bb.chromosome_tree_header.block_size) ||
                bb.btree_leaf_nodes[i].key != bb.key
                break
            end

            return (bb.btree_leaf_nodes[i].chrom_id,
                    bb.btree_leaf_nodes[i].chrom_size)
        end
    end
    error(string("Seqname \"", seqname, "\" is not present in the BigBed file."))
end


function get(bb::BigBedData, seqname::String, first::Int64, last::Int64)
    chrom_id, chrom_size = lookup_seqname(bb, seqname)

    # TODO: Search the R-Tree!
end



