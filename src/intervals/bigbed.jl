# BigBed
# ======
#
# The BigBed format is documented in
#   Kent, W. James, et al. "BigWig and BigBed: # enabling browsing of large
#   distributed datasets." Bioinformatics 26.17 (2010): 2204-2207.
#
# The low level details are documented in a series of tables in the supplement
# of that paper. The immutable types defined below are exactly the data layout
# described in those tables. They are labeled with the table they correspand to.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

const BIGBED_MAGIC = 0x8789F2EB
const BIGWIG_MAGIC = 0x888FFC26
const BIGBED_BTREE_MAGIC = 0x78CA8C91
const BIGBED_RTREE_MAGIC = 0x2468ACE0

const BIGBED_MAX_ZOOM_LEVELS = 10
const BIGBED_CURRENT_VERSION = 4

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

const BIGWIG_DATATYPE_BEDGRAPH  = UInt8(1)
const BIGWIG_DATATYPE_VARSTEP   = UInt8(2)
const BIGWIG_DATATYPE_FIXEDSTEP = UInt8(3)

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

immutable BigBed <: FileFormat end
immutable BigWig <: FileFormat end

type BigBedData <: IntervalStream{BEDMetadata}
    stream::BufferedInputStream
    header::BigBedHeader
    zoom_headers::Vector{BigBedZoomHeader}
    autosql::AbstractString
    summary::BigBedTotalSummary
    btree_header::BigBedBTreeHeader
    rtree_header::BigBedRTreeHeader
    data_count::UInt32

    # preallocated space for reading and searching the B-tree
    btree_internal_nodes::Vector{BigBedBTreeInternalNode}
    btree_leaf_nodes::Vector{BigBedBTreeLeafNode}
    key::Vector{UInt8}
    node_keys::Vector{Vector{UInt8}}
    uncompressed_data::Vector{UInt8}
end

type BigBedDataParser <: AbstractParser
    state::Ragel.State

    # intermediate values used during parsing
    chrom_id::UInt32
    red::Float32
    green::Float32
    blue::Float32
    block_size_idx::Int
    block_first_idx::Int
    seq_names::Nullable{Vector{StringField}}
    assumed_seqname::Nullable{StringField}

    function BigBedDataParser(input::BufferedInputStream;
                              seq_names::Nullable{Vector{StringField}}=Nullable{Vector{StringField}}(),
                              assumed_seqname::Nullable{StringField}=Nullable{StringField}())
        cs = _bigbedparser_start
        return new(Ragel.State(cs, input), 0, 0.0, 0.0, 0.0, 1, 1, seq_names, assumed_seqname)
    end
end

include("bigbed-parser.jl")

"""
Open a BigBed file for reading.

Once opened, entries can be read from the file either by iterating over it, or
by indexing into it with an interval.
"""
function Base.open(stream::BufferedInputStream, ::Type{BigBed})
    # header
    header = read(stream, BigBedHeader)
    if header.magic != BIGBED_MAGIC
        error("Input is not a valid BigBed file.")
    end

    if header.version < 3
        error("Input is an older unsupported version of BigBed.")
    end

    # zoom headers
    zoom_headers = Array(BigBedZoomHeader, header.zoom_levels)
    for i in 1:header.zoom_levels
        zoom_headers[i] = read(stream, BigBedZoomHeader)
    end

    # autosql
    if header.auto_sql_offset == 0
        autosql = ""
    else
        seek(stream, header.auto_sql_offset)
        autosql_buf = IOBuffer()
        while (c = read(stream, UInt8)) != 0x00
            write(autosql_buf, c)
        end
        # TODO: eventually we should parse this and do something useful with it
        autosql = takebuf_string(autosql_buf)
    end

    # total summary
    seek(stream, header.total_summary_offset)
    summary = read(stream, BigBedTotalSummary)

    # b-tree header
    seek(stream, header.chromosome_tree_offset)
    btree_header = read(stream, BigBedBTreeHeader)
    if btree_header.magic != BIGBED_BTREE_MAGIC
        error("BigBed B-Tree magic number was incorrect. File may be corrupt or malformed.")
    end

    # data_count
    seek(stream, header.full_data_offset)
    data_count = read(stream, UInt32)

    # r-tree header
    seek(stream, header.full_index_offset)
    rtree_header = read(stream, BigBedRTreeHeader)
    if rtree_header.magic != BIGBED_RTREE_MAGIC
        error("BigBed R-Tree magic number was incorrect. File may be corrupt or malformed.")
    end

    return BigBedData(stream, header, zoom_headers, autosql, summary,
                      btree_header, rtree_header, data_count,
                      BigBedBTreeInternalNode[BigBedBTreeInternalNode(btree_header.key_size)
                                              for _ in 1:btree_header.block_size],
                      BigBedBTreeLeafNode[BigBedBTreeLeafNode(btree_header.key_size)
                                          for _ in 1:btree_header.block_size],
                      Array(UInt8, btree_header.key_size),
                      Array(Vector{UInt8}, btree_header.block_size),
                      Array(UInt8, header.uncompress_buf_size))
end

function memisless(a::Vector{UInt8}, b::Vector{UInt8})
    if length(a) != length(b)
        return length(a) < length(b)
    end
    i = 1
    while i ≤ length(a)
        if a[i] != b[i]
            return a[i] < b[i]
        end
        i += 1
    end
    return false
end

"""
Return all sequence (name, id, size) tuples in a BigBed B-tree.
"""
function first_btree_leaf_position(bb::BigBedData)
    # find the first leaf-node in the b-tree
    offset = bb.header.chromosome_tree_offset + sizeof(BigBedBTreeHeader)
    seek(bb.stream, offset)
    leafpos = 0
    while true
        node = read(bb.stream, BigBedBTreeNode)
        if node.isleaf != 0
            leafpos = offset
            break
        else
            for i in 1:bb.btree_header.block_size
                read!(bb.stream, bb.btree_internal_nodes[i])
                bb.node_keys[i] = bb.btree_internal_nodes[i].key
            end
            offset = bb.btree_internal_nodes[1].child_offset
            seek(bb.stream, offset)
        end
    end

    if leafpos == 0
        error("Malformed BigBed file: no leaf nodes found in index.")
    end

    return leafpos
end

"""
An iterator over all entries in a BigBed file.
"""
type BigBedIteratorState
    seq_names::Vector{StringField}
    data_count::Int
    data_num::Int
    data_offset::UInt
    parser::BigBedDataParser
    parser_isdone::Bool
    next_interval::Interval{BEDMetadata}
end

function Base.start(bb::BigBedData)
    # read sequence names
    leafpos = first_btree_leaf_position(bb)
    seq_names = Array(StringField, bb.btree_header.item_count)
    leafnode = BigBedBTreeLeafNode(bb.btree_header.key_size)
    seek(bb.stream, leafpos)
    i = 1
    while i <= bb.btree_header.item_count
        node = read(bb.stream, BigBedBTreeNode)
        @assert node.isleaf != 0
        for j in 1:node.count
            read!(bb.stream, leafnode)
            p = findfirst(leafnode.key, 0) - 1
            if p == -1
                p = length(leafnode.key)
            end
            seq_names[i] = StringField(leafnode.key[1:p])
            i += 1
        end
    end

    # read the first data block
    seek(bb.stream, bb.header.full_data_offset)
    data_count = read(bb.stream, UInt64)
    zlib_stream = ZlibInflateInputStream(bb.stream, reset_on_end=false)
    unc_block_size = readbytes!(zlib_stream, bb.uncompressed_data,
                               length(bb.uncompressed_data))

    parser = BigBedDataParser(
        BufferedInputStream(bb.uncompressed_data, unc_block_size),
        seq_names=Nullable(seq_names))

    next_interval = BEDInterval()
    parser_isdone = isnull(tryread!(parser, next_interval))

    return BigBedIteratorState(seq_names, data_count, 1,
                               zlib_stream.source.zstream.total_in,
                               parser, parser_isdone, next_interval)
end

function Base.next(bb::BigBedData, state::BigBedIteratorState)
    value = copy(state.next_interval)

    state.data_num += 1
    if state.data_num > state.data_count
        return value, state
    end

    state.parser_isdone = isnull(tryread!(state.parser, state.next_interval))

    if state.parser_isdone
        seek(bb.stream, bb.header.full_data_offset + state.data_offset + sizeof(UInt64))
        zlib_stream = ZlibInflateInputStream(bb.stream, reset_on_end=false)

        unc_block_size = readbytes!(zlib_stream, bb.uncompressed_data,
                                    length(bb.uncompressed_data))
        state.parser = BigBedDataParser(
             BufferedInputStream(bb.uncompressed_data, unc_block_size),
             seq_names=Nullable(state.seq_names))
        state.data_offset += zlib_stream.source.zstream.total_in

        state.parser_isdone = isnull(tryread!(state.parser, state.next_interval))
        @assert !state.parser_isdone
    end

    return value, state
end

function Base.done(bb::BigBedData, state::BigBedIteratorState)
    return state.data_num > state.data_count
end

"""
An iterator over entries in a BigBed file that intersect a given interval.

Constructed by indexing into `BigBedData` with an interval.
"""
type BigBedIntersectIterator
    bb::BigBedData

    query_seqname::StringField
    query_first::Int64
    query_last::Int64
    query_chrom_id::UInt32
    query_chrom_size::UInt32

    # (offset, size) pairs giving the extents of data blocks that may contain
    # (overlapping intervals
    blocks::Vector{Tuple{UInt64, UInt64}}

    # Index of block currently being parsed
    block_num::Int

    parser::Nullable{BigBedDataParser}
    nextinterval::BEDInterval
    done::Bool
end

"""
Find the given seqname in the BigBed file's index and read the corresponding
sequence id and length
"""
function lookup_seqname(bb::BigBedData, seqname::AbstractString)
    seek(bb.stream, bb.header.chromosome_tree_offset + sizeof(BigBedBTreeHeader))

    fill!(bb.key, 0)
    copy!(bb.key, seqname)

    while true
        node = read(bb.stream, BigBedBTreeNode)
        if node.isleaf == 0
            for i in 1:bb.btree_header.block_size
                read!(bb.stream, bb.btree_internal_nodes[i])
                bb.node_keys[i] = bb.btree_internal_nodes[i].key
            end
            i = searchsortedfirst(bb.node_keys, bb.key, Base.Order.Lt(memisless))
            if !(1 <= i <= bb.btree_header.block_size)
                break
            end

            seek(bb.stream, bb.btree_internal_nodes[i].child_offset)
        else
            for i in 1:bb.btree_header.block_size
                read!(bb.stream, bb.btree_leaf_nodes[i])
                bb.node_keys[i] = bb.btree_leaf_nodes[i].key
            end
            i = searchsortedfirst(bb.node_keys, bb.key, Base.Order.Lt(memisless))
            if !(1 <= i <= bb.btree_header.block_size) ||
                bb.btree_leaf_nodes[i].key != bb.key
                break
            end

            return (bb.btree_leaf_nodes[i].chrom_id,
                    bb.btree_leaf_nodes[i].chrom_size)
        end
    end
    error(string("Seqname \"", seqname, "\" is not present in the BigBed file."))
end

function Base.intersect(bb::BigBedData, query::Interval)
    seqname, first, last = query.seqname, query.first, query.last
    chrom_id, chrom_size = lookup_seqname(bb, seqname)

    # stack of nodes yet to visit
    root_position = bb.header.full_index_offset + sizeof(BigBedRTreeHeader)
    node_positions = UInt64[root_position]

    # stack of data blocks that may contain overlapping intervals
    blocks = Tuple{UInt64, UInt64}[]

    seek(bb.stream, bb.header.full_index_offset + sizeof(BigBedRTreeHeader))
    while !isempty(node_positions)
        pos = pop!(node_positions)
        seek(bb.stream, pos)
        node = read(bb.stream, BigBedRTreeNode)
        if node.isleaf == 0
            for i in 1:node.count
                internal_node = read(bb.stream, BigBedRTreeInternalNode)
                if internal_node.start_chrom_ix <= chrom_id <= internal_node.end_chrom_ix &&
                    (chrom_id < internal_node.end_chrom_ix || first <= internal_node.end_base) &&
                    (chrom_id > internal_node.start_chrom_ix || last - 1 >= internal_node.start_base)
                    push!(node_positions, internal_node.data_offset)
                end
            end
        else
            for i in 1:node.count
                leaf_node = read(bb.stream, BigBedRTreeLeafNode)
                if leaf_node.start_chrom_ix <= chrom_id <= leaf_node.end_chrom_ix &&
                    (chrom_id < leaf_node.end_chrom_ix || first <= leaf_node.end_base) &&
                    (chrom_id > leaf_node.start_chrom_ix || last - 1 >= leaf_node.start_base)
                    push!(blocks, (leaf_node.data_offset, leaf_node.data_size))
                end
            end
        end
    end

    return BigBedIntersectIterator(bb, seqname, first, last, chrom_id,
                                   chrom_size, blocks, 0,
                                   Nullable{BigBedDataParser}(),
                                   BEDInterval(), false)
end

function find_next_intersection!(it::BigBedIntersectIterator)
    it.done = true
    while !isnull(it.parser) || it.block_num < length(it.blocks)
        it.done = true
        while !isnull(it.parser)
            parser = get(it.parser)
            # Run the parser until we find an intersection
            it.done = isnull(tryread!(parser, it.nextinterval))
            if it.done
                it.parser = Nullable{BigBedDataParser}()
                break
            end

            if parser.chrom_id == it.query_chrom_id &&
               it.nextinterval.first <= it.query_last &&
               it.nextinterval.last >= it.query_first
                it.done = false
                return
            end
        end

        if it.block_num < length(it.blocks)
            # advance to the next block of interest ant initialize a new parser
            it.block_num += 1
            block_offset, block_size = it.blocks[it.block_num]

            seek(it.bb.stream, block_offset)
            @assert block_size <= length(it.bb.uncompressed_data)
            unc_block_size = readbytes!(ZlibInflateInputStream(it.bb.stream, reset_on_end=false),
                                        it.bb.uncompressed_data,
                                        length(it.bb.uncompressed_data))
            it.parser = BigBedDataParser(
                BufferedInputStream(it.bb.uncompressed_data, unc_block_size),
                assumed_seqname=Nullable(it.query_seqname))
        else
            it.done = true
            return
        end
    end
end

function Base.start(it::BigBedIntersectIterator)
    find_next_intersection!(it)
    return nothing
end

function Base.next(it::BigBedIntersectIterator, ::Void)
    value = copy(it.nextinterval)
    find_next_intersection!(it)
    return value, nothing
end

function Base.done(it::BigBedIntersectIterator, ::Void)
    return it.done
end


# BigBed Output
# -------------

function write_zeros(out::IO, n::Integer)
    for i in 1:n
        write(out, 0x00)
    end
end

immutable BigBedChromInfo
    name::StringField   # chromosome name
    id::UInt32          # unique id for chromosome
    size::UInt32        # size of chromosome
    item_count::UInt64  # number of items
end

immutable BigBedBounds
    offset::UInt64
    chrom_ix::UInt32
    start::UInt32
    stop::UInt32
end

function bigbed_chrom_info(intervals::IntervalCollection, chrom_sizes::Dict)
    info = Array(BigBedChromInfo, length(intervals.trees))
    seqnames = collect(StringField, keys(intervals.trees))
    sort!(seqnames)
    for (i, seqname) in enumerate(seqnames)
        info[i] = BigBedChromInfo(
            seqname,
            i - 1,
            chrom_sizes[seqname],
            length(intervals.trees[seqname]))
    end
    return info
end

function bigbed_btree_count_levels(max_block_size, item_count)
    levels = 1
    while item_count > max_block_size
        item_count = div(item_count + max_block_size - 1, max_block_size)
        levels += 1
    end
    return levels
end

function bigbed_btree_write_index_level(out::IO, block_size,
                                        items::Vector{BigBedChromInfo},
                                        key_size, index_offset, level)
    # See writeIndexLevel in bPlusTree.c

    # calculate number of nodes to write at this level
    slot_size_per = block_size ^ level
    node_size_per = slot_size_per * block_size
    node_count = div(length(items) + node_size_per - 1, node_size_per)

    # calculate sizes and offsets
    block_header_size = 4
    value_size = 8
    bytes_in_index_block = block_header_size + block_size * (key_size + sizeof(UInt64))
    bytes_in_leaf_block = block_header_size + block_size * (key_size + value_size)
    bytes_in_next_level_block = level == 1 ? bytes_in_leaf_block : bytes_in_index_block
    level_size = node_count * bytes_in_index_block
    end_level = index_offset + level_size
    next_child = end_level

    for i in 1:node_size_per:length(items)
        # calculate size of this block
        count_one = min(block_size, div(length(items) - (i-1) + slot_size_per - 1, slot_size_per))

        # write block header
        write(out, BigBedBTreeNode(0, 0, count_one))

        # write used slots
        slots_used = 0
        end_ix = min(length(items), i + node_size_per)
        for j in i:slot_size_per:end_ix
            item = items[j]
            write(out, item.name)
            # pad name
            for k in length(item.name)+1:key_size
                write(out, 0x00)
            end
            write(out, convert(UInt64, next_child))
            next_child += bytes_in_next_level_block
            slots_used += 1
        end
        @assert slots_used == count_one "Incorrect number of items written to index"

        # write out empty slots as all zero
        slot_size = key_size + 8
        for j in count_one:(block_size-1)
            for k in 1:slot_size
                write(out, 0x00)
            end
        end
    end

    return end_level
end

function bigbed_btree_write_leaf_level(out::IO, block_size,
                                       items::Vector{BigBedChromInfo},
                                       key_size, val_size)
    # See writeLeafLevel in bPlusTree.c

    count_one = 0
    count_left = length(items)
    i = 0
    while i < length(items)
        count_one = min(count_left, block_size)

        # write block header
        write(out, BigBedBTreeNode(1, 0, count_one))

        # write out position in genome and in file for each item
        for j in 1:count_one
            @assert i + j <= length(items) "Incorrect item count when writing B-tree leaves"
            item = items[i + j]
            write(out, item.name)
            # pad name
            for k in length(item.name)+1:key_size
                write(out, 0x00)
            end
            write(out, convert(UInt32, item.id))
            write(out, convert(UInt32, item.size))
        end

        # pad out any unused bits of last block with zeroes
        slot_size = key_size + val_size
        for j in count_one:(block_size-1)
            for k in 1:slot_size
                write(out, 0x00)
            end
        end

        count_left -= count_one
        i += count_one
    end
end

function bigbed_write_chrom_tree(out::IO, intervals::IntervalCollection,
                                 block_size,
                                 chrom_sizes::Dict=Dict{AbstractString, Int64}())
    # See bbiWriteChromInfo in bbiWrite.c
    length(intervals.trees)
    for seqname in keys(intervals.trees)
        if !haskey(chrom_sizes, seqname)
            chrom_sizes[seqname] =
                length(intervals.trees[seqname]) > 0 ?
                    intervals.trees[seqname].root.maxend : 0
        end
    end

    max_chrom_name_size = 0
    for seqname in keys(intervals.trees)
        max_chrom_name_size = max(max_chrom_name_size, length(seqname))
    end

    info = bigbed_chrom_info(intervals, chrom_sizes)
    item_count = length(info)
    chrom_block_size = min(block_size, item_count)
    value_size = 8; # sizeof(chrom_id) + sizeof(chrom_size)

    # write header. See bptFileBulkIndex in bPlusTree.c
    header = BigBedBTreeHeader(
        BIGBED_BTREE_MAGIC,
        chrom_block_size,
        max_chrom_name_size,
        value_size,
        item_count,
        0)
    write(out, header)
    index_offset = position(out)

    # write non-leaf nodes.
    levels = bigbed_btree_count_levels(chrom_block_size, item_count)
    for i in (levels-1):-1:1
        end_level_offset = bigbed_btree_write_index_level(out, chrom_block_size,
                                                          info, max_chrom_name_size,
                                                          index_offset, i)
        index_offset = position(out)
        @assert index_offset == end_level_offset "Incorrect file offsets reported."
    end

    # write leaf nodes
    bigbed_btree_write_leaf_level(out, chrom_block_size, info,
                                  max_chrom_name_size, value_size)
    return info
end

function bigbed_count_sections_needed(items::Vector{BigBedChromInfo}, items_per_slot)
    count = 0
    for item in items
        count_one = div(item.item_count + items_per_slot - 1, items_per_slot)
        count += count_one
    end
    return count
end

function bigbed_write_blocks(out::IO, intervals::IntervalCollection,
                             chrom_info::Vector{BigBedChromInfo},
                             items_per_slot, bounds::Vector{BigBedBounds},
                             section_count, compressed::Bool,
                             res_try_count, res_scales, res_sizes)
    # See writeBlocks in bedToBigBed.c

    section_ix = 0
    buf = IOBuffer()
    res_ends = zeros(Int, res_try_count)
    max_block_size = 0

    for info in chrom_info
        tree = intervals.trees[info.name]
        it_state = start(tree)
        fill!(res_ends, 0)

        while !done(tree, it_state)
            section_ix += 1
            item_ix = 1
            block_offset = position(out)
            first_interval = first(intervals.trees[info.name])
            start_pos = -1
            end_pos  = first_interval.last

            while item_ix <= items_per_slot && !done(tree, it_state)
                interval, it_state = next(tree, it_state)

                interval_start = interval.first - 1 # bed is 0-based
                interval_end   = interval.last
                if start_pos < 0
                    start_pos = interval_start
                end
                end_pos = max(end_pos, interval_end)

                write(buf, convert(UInt32, info.id))
                write(buf, convert(UInt32, interval_start))
                write(buf, convert(UInt32, interval_end))
                write_optional_fields(buf, interval, false)
                write(buf, '\0')

                item_ix += 1

                for res_try in 1:res_try_count
                    res_end = res_ends[res_try]
                    if interval_start >= res_end
                        res_sizes[res_try] += 1
                        res_ends[res_try] = res_end = interval_start + res_scales[res_try]
                    end
                    while interval_end > res_end
                        res_sizes[res_try] += 1
                        res_ends[res_try] = res_end = res_end + res_scales[res_try]
                    end
                end
            end

            block = takebuf_array(buf)
            max_block_size = max(max_block_size, length(block))

            if compressed
                writer = ZlibDeflateOutputStream(out)
                write(writer, block)
                flush(writer)
            else
                write(out, block)
            end

            bounds[section_ix] = BigBedBounds(block_offset, info.id,
                                              start_pos, end_pos)
        end
    end

    @assert section_ix == section_count "Incorrect number of sections when writing BigBed data"

    return max_block_size
end

function bigwig_write_blocks{T<:Number}(out::IO, intervals::IntervalCollection{T},
                                        chrom_info::Vector{BigBedChromInfo},
                                        items_per_slot, bounds::Vector{BigBedBounds},
                                        section_count, compressed::Bool,
                                        res_try_count, res_scales, res_sizes)
    section_ix = 0
    res_ends = zeros(Int, res_try_count)
    items = Array(BedGraphItem, items_per_slot)
    max_section_size = 0

    # See writeSections in BedGraphToBigWig.c
    for info in chrom_info
        tree = intervals.trees[info.name]
        it_state = start(tree)
        fill!(res_ends, 0)

        while !done(tree, it_state)
            section_ix += 1
            item_ix = 1

            while item_ix <= items_per_slot && !done(tree, it_state)
                interval, it_state = next(tree, it_state)
                interval_start = interval.first - 1
                interval_end = interval.last
                items[item_ix] = BedGraphItem(
                    interval_start, interval_end, interval.metadata)
                item_ix += 1

                for res_try in 1:res_try_count
                    res_end = res_ends[res_try]
                    if interval_start >= res_end
                        res_sizes[res_try] += 1
                        res_ends[res_try] = res_end = interval_start + res_scales[res_try]
                    end
                    while interval_end > res_end
                        res_sizes[res_try] += 1
                        res_ends[res_try] = res_end = res_end + res_scales[res_try]
                    end
                end
            end

            # write a section
            chrom_id = info.id
            section_start = items[1].chrom_start
            section_end = items[end].chrom_end

            section_offset = position(out)
            header = BigWigSectionHeader(chrom_id, section_start, section_end,
                                         0, 0, BIGWIG_DATATYPE_BEDGRAPH, 0, item_ix - 1)

            if compressed
                writer = ZlibDeflateOutputStream(out)
                write(writer, header)
                for i in 1:item_ix-1
                    write(writer, items[i])
                end
                flush(writer)
            else
                write(out, header)
                for i in 1:item_ix-1
                    write(out, items[i])
                end
            end
            section_end_offset = position(out)
            max_section_size = max(max_section_size, section_end_offset - section_offset)

            bounds[section_ix] = BigBedBounds(section_offset, info.id,
                                              section_start, section_end)
        end
    end

    @assert section_ix == section_count "Incorrect number of sections when writing BigWig data"

    return max_section_size
end

type BigBedRTree
    next::Nullable{BigBedRTree}
    children::Nullable{BigBedRTree}
    parent::Nullable{BigBedRTree}
    start_chrom_ix::UInt32
    start_base::UInt32
    end_chrom_ix::UInt32
    end_base::UInt32
    start_file_offset::UInt64
    end_file_offset::UInt64

    function BigBedRTree(next, children, parent, start_chrom_ix, start_base,
                         end_chrom_ix, end_base, start_file_offset, end_file_offset)
        return new(next, children, parent, start_chrom_ix, end_chrom_ix,
                   start_base, end_base, start_file_offset, end_file_offset)
    end

    function BigBedRTree()
        return new()
    end
end

function Base.copy(tree::BigBedRTree)
    return BigBedRTree(tree.next, tree.children, tree.parent,
                       tree.start_chrom_ix, tree.start_base,
                       tree.end_chrom_ix, tree.end_base,
                       tree.start_file_offset, tree.end_file_offset)
end

function Base.length(list::Nullable{BigBedRTree})
    count = 0
    while !isnull(list)
        count += 1
        list = get(list).next
    end
    return count
end

function bigbed_calc_rtree_level_sizes!(level_sizes, tree::Nullable{BigBedRTree},
                                        level, max_level)
    # See calcLevelSizes in cirTree.c

    el = tree
    while !isnull(el)
        level_sizes[level] += 1
        if level <= max_level
            bigbed_calc_rtree_level_sizes!(level_sizes, get(el).children,
                                           level + 1, max_level)
        end
        el = get(el).next
    end
end

function bigbed_rtree_write_index_level(out::IO, block_size, child_node_size,
                                        tree::BigBedRTree, offset_of_first_child,
                                        cur_level, dest_level)
    offset = offset_of_first_child
    if cur_level == dest_level
        # we've reached the right level, write out a node header
        count_one = length(tree.children)
        write(out, BigBedRTreeNode(0, 0, count_one))

        # write out elements of this node
        nel = tree.children
        while !isnull(nel)
            el = get(nel)
            write(out, BigBedRTreeInternalNode(el.start_chrom_ix, el.start_base,
                                               el.end_chrom_ix, el.end_base, offset))
            offset += child_node_size
            nel = el.next
        end

        # write out zeroes for empty slots in node
        for i in count_one:(block_size-1)
            for j in 1:sizeof(BigBedRTreeInternalNode)
                write(out, 0x00)
            end
        end
    else
        # otherwise recurse on children
        nel = tree.children
        while !isnull(nel)
            el = get(nel)
            offset = bigbed_rtree_write_index_level(out, block_size,
                                                    child_node_size, el,
                                                    cur_level + 1,
                                                    dest_level,
                                                    offset)
            nel = el.next
        end
    end
    return offset
end

function bigbed_rtree_write_leaf_level(out::IO, items_per_slot, leaf_node_size,
                                       tree::BigBedRTree, cur_level, leaf_level)
    if cur_level == leaf_level
        # we've reached the right level, write out node header
        count_one = length(tree.children)
        write(out, BigBedRTreeNode(1, 0, count_one))

        nel = tree.children
        while !isnull(nel)
            el = get(nel)
            write(out, BigBedRTreeLeafNode(el.start_chrom_ix, el.start_base,
                                           el.end_chrom_ix, el.end_base,
                                           el.start_file_offset,
                                           el.end_file_offset - el.start_file_offset))
            nel = el.next
        end

        # write out zeroes for empty slots in node
        for i in count_one:(items_per_slot-1)
            for j in 1:sizeof(BigBedRTreeInternalNode)
                write(out, 0x00)
            end
        end
    else
        # otherwise recurse on children
        nel = tree.children
        while !isnull(nel)
            el = get(nel)
            bigbed_rtree_write_leaf_level(out, items_per_slot, leaf_node_size,
                                          el, cur_level + 1, leaf_level)
            nel = el.next
        end
    end
end

function Base.write(out::IO, tree::BigBedRTree, block_size, level_count)
    # See writeTreeToOpenFile in cirTree.c

    level_sizes = zeros(Int, level_count)
    bigbed_calc_rtree_level_sizes!(level_sizes, Nullable{BigBedRTree}(tree), 1, level_count)

    level_offsets = Array(UInt64, level_count)

    index_node_size = sizeof(BigBedRTreeNode) + block_size * sizeof(BigBedRTreeInternalNode)
    leaf_node_size = sizeof(BigBedRTreeNode) + block_size * sizeof(BigBedRTreeLeafNode)

    offset = position(out)
    for i in 1:level_count
        level_offsets[i] = offset
        offset += level_sizes[i] * index_node_size
    end

    final_level = level_count - 2
    for i in 1:final_level
        child_node_size = i == final_level ? leaf_node_size : index_node_size
        bigbed_rtree_write_index_level(out, block_size, child_node_size, tree,
                                       level_offsets[i+1], 1, i)
        @assert position(out) == level_offsets[i+1]
    end

    leaf_level = level_count - 1
    bigbed_rtree_write_leaf_level(out, block_size, leaf_node_size, tree,
                                  1, leaf_level)
end

function Base.reverse!(tree::Nullable{BigBedRTree})
    if isnull(tree)
        return tree
    end

    nodes = BigBedRTree[]
    node = tree
    while !isnull(node)
        push!(nodes, get(node))
        node = get(node).next
    end

    for i in length(nodes):-1:2
        nodes[i].next = nodes[i-1]
    end
    nodes[1].next = Nullable{BigBedRTree}()
    return Nullable{BigBedRTree}(nodes[end])
end

function bigbed_build_rtree(bounds::Vector{BigBedBounds}, block_size,
                            items_per_slot, end_file_offset)
    # See rTreeFromChromRangeArray in cirTree.c

    list  = Nullable{BigBedRTree}()
    next_offset = isempty(bounds) ? 0 : bounds[1].offset
    for i in 1:items_per_slot:length(bounds)
        final_iteration = false
        one_size = length(bounds) - i + 1
        if one_size > items_per_slot
            one_size = items_per_slot
        else
            final_iteration = true
        end

        bound = bounds[i]
        el = BigBedRTree()
        el.next = Nullable{BigBedRTree}()
        el.children = Nullable{BigBedRTree}()
        el.parent = Nullable{BigBedRTree}()
        el.start_chrom_ix = el.end_chrom_ix = bound.chrom_ix
        el.start_base = bound.start
        el.end_base = bound.stop
        el.start_file_offset = next_offset

        # Figure out end of element from offset of next element (or file size
        # for final element.)
        if final_iteration
            next_offset = end_file_offset
        else
            next_offset = bounds[i + one_size].offset
        end
        el.end_file_offset = next_offset

        # expand area spanned to include all items in block
        for j in 2:one_size
            bound = bounds[i + j - 1]
            if bound.chrom_ix < el.start_chrom_ix
                el.start_chrom_ix = bound.chrom_ix
                el.start_base = bound.start
            elseif bound.chrom_ix == el.start_chrom_ix
                el.start_base = min(el.start_base, bound.start)
            elseif bound.chrom_ix > el.end_chrom_ix
                el.end_chrom_ix = bound.chrom_ix
                el.end_base = bound.stop
            elseif bound.chrom_ix == el.end_chrom_ix
                el.end_base = max(el.end_base, bound.stop)
            end
        end

        el.next = list
        list = Nullable{BigBedRTree}(el)
    end
    list = reverse!(list)

    # Now iterate through making more and more condensed versions until have
    # just one
    tree = list
    level_count = 1
    while (!isnull(tree) && !isnull(get(tree).next)) || level_count < 2
        list = Nullable{BigBedRTree}()
        slots_used = block_size
        parent = Nullable{BigBedRTree}()
        nel = tree
        while !isnull(nel)
            el = get(nel)
            next_el = el.next
            if slots_used >= block_size
                slots_used = 1
                parent = Nullable{BigBedRTree}(copy(el))
                get(parent).children = el
                el.parent = parent
                el.next = Nullable{BigBedRTree}()

                get(parent).next = list
                list = parent
            else
                slots_used += 1
                el.next = get(parent).children
                get(parent).children = el
                el.parent = parent
                if el.start_chrom_ix < get(parent).start_chrom_ix
                    get(parent).start_chrom_ix = el.start_chrom_ix
                    get(parent).start_base = el.start_base
                elseif el.start_chrom_ix == get(parent).start_chrom_ix
                    get(parent).start_base = min(el.start_base, get(parent).start_base)
                end

                if el.end_chrom_ix > get(parent).end_chrom_ix
                    get(parent).end_chrom_ix = el.end_chrom_ix
                    get(parent).end_base = el.end_base
                elseif el.end_chrom_ix == get(parent).end_chrom_ix
                    get(parent).end_base  = max(el.end_base, get(parent).end_base)
                end
            end
            nel = next_el
        end

        list = reverse!(list)
        el = list
        while !isnull(el)
            get(el).children = reverse!(get(el).children)
            el = get(el).next
        end
        tree = list
        level_count += 1
    end

    @assert !isnull(tree)

    (get(tree), level_count)
end

function bigbed_write_index(out::IO, bounds::Vector{BigBedBounds},
                            block_size, items_per_slot, end_file_offset)
    # See cirTreeFileBulkIndexToOpenFile in cirTree.c

    rtree, level_count = bigbed_build_rtree(bounds, block_size, items_per_slot,
                                            end_file_offset)
    write(out, BigBedRTreeHeader(BIGBED_RTREE_MAGIC, block_size, length(bounds),
                                 rtree.start_chrom_ix, rtree.start_base,
                                 rtree.end_chrom_ix, rtree.end_base,
                                 end_file_offset, items_per_slot, 0))
    write(out, rtree, block_size, level_count)
end

# Write out sum to file, keeping track of minimal info on it in, and also adding
# it to second level summary
function output_one_summary_further_reduce(out::IO, sum::BigBedZoomData,
                                           twice_reduced_list::Vector{BigBedZoomData},
                                           double_reduction_size,
                                           bounds_array::Vector{BigBedBounds},
                                           bounds_array_idx,
                                           chrom_size)
    # See bbiOutputOneSummaryFurtherReduce in bbiWrite.c

    @assert bounds_array_idx <= length(bounds_array)

    bounds_array[bounds_array_idx] = BigBedBounds(position(out), sum.chrom_id,
                                                  sum.chrom_start, sum.chrom_end)
    write(out, sum)

    # Fold summary info into twice_reduced_list
    if isempty(twice_reduced_list) || twice_reduced_list[end].chrom_id != sum.chrom_id ||
       twice_reduced_list[end].chrom_start + double_reduction_size < sum.chrom_end
        push!(twice_reduced_list, sum)
    else
        twice_reduced_list[end] = BigBedZoomData(
            twice_reduced_list[end].chrom_id,
            twice_reduced_list[end].chrom_start,
            sum.chrom_end,
            twice_reduced_list[end].valid_count + sum.valid_count,
            min(sum.min_val, twice_reduced_list[end].min_val),
            max(sum.max_val, twice_reduced_list[end].max_val),
            twice_reduced_list[end].sum_data + sum.sum_data,
            twice_reduced_list[end].sum_squares + sum.sum_squares)
    end

    return bounds_array_idx + 1
end

# Write out data reduced by factor of initialReduction.  Also calculate and keep
# in memory next reduction level.  This is more work than some ways, but it
# keeps us from having to keep the first reduction entirely in memory.
function bigbed_write_reduced_once_return_reduced_twice(
                out::IO, intervals::IntervalCollection,
                interval_tree_transform::Function,
                chrom_info::Vector{BigBedChromInfo},
                initial_reduction, initial_reduced_count,
                zoom_increment, block_size, items_per_slot, compressed)

    # See writeReducedOnceReturnReducedTwice in bedToBigBed.c

    double_reduction_size = initial_reduction * zoom_increment
    bounds_array = Array(BigBedBounds, initial_reduced_count)
    bounds_array_idx = 1
    twice_reduced_list = BigBedZoomData[]

    ret_data_start = position(out)
    write(out, convert(UInt32, initial_reduced_count))

    # This gets a little complicated I'm afraid.  The strategy is to:
    #  1) Build up a range tree that represents coverage depth on that chromosome
    #     This also has the nice side effect of getting rid of overlaps.
    #  2) Stream through the range tree, outputting the initial summary level and
    #     further reducing.
    #

    valid_count = UInt64(0)
    min_val = 0.0
    max_val = 0.0
    sum_data = 0.0
    sum_squares = 0.0
    first_time = true

    sum = Nullable{BigBedZoomData}()

    for info in chrom_info
        sum = Nullable{BigBedZoomData}()
        for range in interval_tree_transform(intervals.trees[info.name])
            val = convert(Float64, range.metadata)
            chrom_start = range.first - 1
            chrom_end = range.last
            size = chrom_end - chrom_start
            @assert size >= 0 "Bad size: $chrom_start, $chrom_end"
            if first_time
                valid_count = size
                min_val = max_val = val
                sum_data = val * size
                sum_squares = val * val * size
                first_time = false
            else
                valid_count += size
                min_val = min(min_val, val)
                max_val = max(max_val, val)
                sum_data += val * size
                sum_squares += val * val * size
            end

            # If start past existing block then output it
            if !isnull(sum) && get(sum).chrom_end <= chrom_start
                bounds_array_idx = output_one_summary_further_reduce(
                        out, get(sum), twice_reduced_list,
                        double_reduction_size, bounds_array,
                        bounds_array_idx, info.size)
                sum = Nullable{BigBedZoomData}()
            end

            # If don't have a summary we're working on now, make one
            if isnull(sum)
                sum = Nullable{BigBedZoomData}(
                    BigBedZoomData(info.id, chrom_start,
                                   min(info.size, chrom_start + initial_reduction),
                                   0, val, val, 0.0, 0.0))
            end

            # Deal with case where might have to split an item between multiple
            # summaries. This loop handles all but the final affected summary in
            # that case
            s = get(sum)
            while chrom_end > s.chrom_end
                # Fold in bits that ovelap with existing summary and output
                overlap = min(chrom_end, s.chrom_end) - max(chrom_start, s.chrom_start)
                @assert overlap > 0
                @assert overlap <= size

                s = BigBedZoomData(
                        s.chrom_id, s.chrom_start, s.chrom_end,
                        s.valid_count + overlap,
                        min(s.min_val, val), max(s.max_val, val),
                        s.sum_data + val * overlap,
                        s.sum_squares + val * val * overlap)

                bounds_array_idx = output_one_summary_further_reduce(
                        out, s, twice_reduced_list, double_reduction_size,
                        bounds_array, bounds_array_idx, info.size)
                size -= overlap

                # Move summary to next part
                chrom_start = s.chrom_end
                #chrom_end = min(info.size, chrom_start + initial_reduction)
                s = BigBedZoomData(
                        s.chrom_id, chrom_start,
                        min(info.size, chrom_start + initial_reduction),
                        0, val, val, 0.0, 0.0)
            end

            # Add to summary
            s = BigBedZoomData(
                    s.chrom_id, s.chrom_start, s.chrom_end,
                    s.valid_count + size,
                    min(s.min_val, val), max(s.max_val, val),
                    s.sum_data + val * size,
                    s.sum_data + val * val * size)
            sum = Nullable(s)
        end

        if !isnull(sum)
            bounds_array_idx = output_one_summary_further_reduce(
                    out, get(sum), twice_reduced_list, double_reduction_size,
                    bounds_array, bounds_array_idx, info.size)
        end
    end

    # Write out 1st zoom index
    index_offset = position(out)
    @assert bounds_array_idx == length(bounds_array) + 1
    bigbed_write_index(out, bounds_array, block_size, items_per_slot,
                       index_offset)
    return twice_reduced_list, BigBedTotalSummary(valid_count, min_val, max_val,
                                                  sum_data, sum_squares)
end

function bigbed_write_summary_and_index(
        out::IO, summary_list::Vector{BigBedZoomData}, block_size,
        items_per_slot, compressed)
    # See: bbiWriteSummaryAndIndex in bbiWrite.c
    if compressed
        bigbed_write_summary_and_index_comp(out, summary_list, block_size, items_per_slot)
    else
        bigbed_write_summary_and_index_unc(out, summary_list, block_size, items_per_slot)
    end
end

function bigbed_write_summary_and_index_unc(
        out::IO, summary_list::Vector{BigBedZoomData}, block_size,
        items_per_slot)
    # See: bbiWriteSummaryAndIndexUnc in bbiWrite.c
    count = length(summary_list)
    write(out, convert(UInt32, count))
    bounds_array = Array(BigBedBounds, count)
    for (i, summary) in enumerate(summary_list)
        bounds_array[i] = BigBedBounds(
            position(out), summary.chrom_id, summary.chrom_start,
            summary.chrom_end)
        write(out, summary)
    end
    index_offset = position(out)
    bigbed_write_index(out, bounds_array, block_size, items_per_slot, false)
    return index_offset
end

function bigbed_write_summary_and_index_comp(
        out::IO, summary_list::Vector{BigBedZoomData}, block_size,
        items_per_slot)
    # See: bbiWriteSummaryAndIndexComp in bbiWrite.c
    count = length(summary_list)
    write(out, convert(UInt32, count))
    bounds_array = Array(BigBedBounds, count)

    items_left = count
    sum_ix = 1
    i = 1
    while items_left > 0
        items_in_slot = min(items_per_slot, items_left)
        writer = ZlibDeflateOutputStream(out)
        file_pos = position(out)
        for _ in 1:items_in_slot
            summary = summary_list[sum_ix]
            bounds_array[sum_ix] = BigBedBounds(file_pos, summary.chrom_id,
                                                summary.chrom_start, summary.chrom_end)
            write(writer, summary_list[sum_ix])
            sum_ix += 1
            if sum_ix > count
                break
            end
        end
        items_left -= items_in_slot
        flush(writer)
    end
    index_offset = position(out)
    bigbed_write_index(out, bounds_array, block_size, items_per_slot, false)
    return index_offset
end

function bigbed_summary_simple_reduce(summary_list::Vector{BigBedZoomData}, reduction)
    new_summary_list = BigBedZoomData[]
    for summary in summary_list
        if isempty(new_summary_list) ||
           new_summary_list[end].chrom_id != summary.chrom_id ||
           summary.chrom_end > new_summary_list[end].chrom_start + reduction
            push!(new_summary_list, summary)
        else
            @assert new_summary_list[end].chrom_end < summary.chrom_end
            new_summary_list[end] = BigBedZoomData(
                new_summary_list[end].chrom_id,
                new_summary_list[end].chrom_start,
                summary.chrom_end,
                new_summary_list[end].valid_count + summary.valid_count,
                min(new_summary_list[end].min_val, summary.min_val),
                max(new_summary_list[end].max_val, summary.max_val),
                new_summary_list[end].sum_data + summary.sum_data,
                new_summary_list[end].sum_squares + summary.sum_squares)
       end
    end
    return new_summary_list
end

function bed_field_count(intervals::IntervalCollection{BEDMetadata})
    num_fields = 0
    for interval in intervals
        interval_num_fields = 3 + interval.metadata.used_fields
        if num_fields == 0 || num_fields > interval_num_fields
            num_fields = interval_num_fields
        end
    end
    return num_fields
end

function Base.write(out::IO, ::Type{BigBed}, intervals::IntervalCollection{BEDMetadata};
                    block_size::Int=256, items_per_slot::Int=512,
                    compressed::Bool=true)
    write_bigbed_bigwig(out, BigBed, intervals, block_size, items_per_slot, compressed)
end

function Base.write{T<:Number}(out::IO, ::Type{BigWig}, intervals::IntervalCollection{T};
                               block_size::Int=256, items_per_slot::Int=512,
                               compressed::Bool=true)
    write_bigbed_bigwig(out, BigWig, intervals, block_size, items_per_slot, compressed)
end

# Unified write function for bigwig and bigbed
function write_bigbed_bigwig(out::IO, fmt::Union{Type{BigBed}, Type{BigWig}},
                             intervals::IntervalCollection,
                             block_size::Int=256, items_per_slot::Int=512,
                             compressed::Bool=true)

    # See the function bbFileCreate in bedToBigBed.c in kent to see the only
    # other implementation of this I'm aware of. Also bedGraphToBigWig in
    # bedGraphToBigWig.c

    if !applicable(seek, out, 1)
        error("BigBed can only be written to seekable output streams.")
    end

    field_count = fmt == BigBed ? bed_field_count(intervals) : 0

    # write dummy headers, come back and fill them later
    write_zeros(out, sizeof(BigBedHeader))
    write_zeros(out, BIGBED_MAX_ZOOM_LEVELS * sizeof(BigBedZoomHeader))

    # TODO: optional autoSql specification
    as_offset = 0

    total_summary_offset = position(out)
    total_summary = BigBedTotalSummary(0, 0.0, 0.0, 0.0, 0.0)
    write(out, total_summary)

    chrom_tree_offset = position(out)
    chrom_info = bigbed_write_chrom_tree(out, intervals, block_size)

    # set up to keep track of possible initial reduction levels
    total_bases = 0
    for interval in intervals
        total_bases += interval.last - interval.first + 1
    end
    ave_span = total_bases / length(intervals)

    res_try_count = 10
    res_increment = 4
    res_scales = Array(Int, res_try_count)
    res_sizes = Array(Int, res_try_count)
    min_zoom = 10
    res = max(round(Int, ave_span), min_zoom)
    for res_try in 1:res_try_count
        res_sizes[res_try] = 0
        res_scales[res_try] = res
        if res > 1000000000
            res_try_count = res_try
            break
        end
        res *= res_increment
    end

    # write out primary full resolution data in sections, collect stats to
    # use for reductions
    data_offset = position(out)
    # TODO: The BigBed/BigWig specification differs from Kent's impementation
    # here. The spec says this count is a 32bit int, but the implmentation uses
    # a 64bit int. We follow the implementation for practical reasons.
    write(out, convert(UInt64, length(intervals)))
    block_count = bigbed_count_sections_needed(chrom_info, items_per_slot)
    bounds = Array(BigBedBounds, block_count)

    if fmt == BigBed
        max_block_size = bigbed_write_blocks(out, intervals, chrom_info, items_per_slot,
                                             bounds, block_count, compressed,
                                             res_try_count, res_scales, res_sizes)
    else
        max_block_size = bigwig_write_blocks(out, intervals, chrom_info, items_per_slot,
                                             bounds, block_count, compressed,
                                             res_try_count, res_scales, res_sizes)
    end

    # write out primary data index
    index_offset = position(out)
    bigbed_write_index(out, bounds, block_size, 1, index_offset)

    # declare arrays and vars that track the zoom levels we actually output
    zoom_levels = 0
    zoom_amounts = Array(UInt32, BIGBED_MAX_ZOOM_LEVELS)
    zoom_data_offsets = Array(UInt64, BIGBED_MAX_ZOOM_LEVELS)
    zoom_index_offsets = Array(UInt64, BIGBED_MAX_ZOOM_LEVELS)
    zoom_levels = 0

    if ave_span > 0
        data_size = index_offset - data_offset
        max_reduced_size = div(data_size, 2)
        initial_reduction = 0
        initial_reduced_count = 0

        for res_try in 1:res_try_count
            reduced_size = res_sizes[res_try] * sizeof(BigBedZoomData)
            if compressed
                reduced_size = div(reduced_size, 2) # estimate!
            end

            if reduced_size <= max_reduced_size
                initial_reduction = res_scales[res_try]
                initial_reduced_count = res_sizes[res_try]
                break
            end
        end

        if initial_reduction > 0
            zoom_increment = 4
            rezoomed_list, total_summary = bigbed_write_reduced_once_return_reduced_twice(
                    out, intervals,
                    fmt == BigBed ? coverage : identity,
                    chrom_info, initial_reduction,
                    initial_reduced_count, zoom_increment, block_size,
                    items_per_slot, compressed)

            zoom_amounts[1] = initial_reduction
            zoom_levels = 1

            zoom_count = initial_reduced_count
            reduction = initial_reduction * zoom_increment
            while zoom_levels < BIGBED_MAX_ZOOM_LEVELS
                rezoom_count = length(rezoomed_list)
                if rezoom_count >= zoom_count
                    break
                end

                zoom_count = rezoom_count
                zoom_data_offsets[zoom_levels] = position(out)
                zoom_index_offsets[zoom_levels] =
                        bigbed_write_summary_and_index(out, rezoomed_list, block_size,
                                                       items_per_slot, compressed)
                zoom_amounts[zoom_levels] = reduction
                zoom_levels += 1
                reduction *= zoom_increment

                rezoomed_list = bigbed_summary_simple_reduce(rezoomed_list, reduction)
            end
        end
    end

    uncompress_buf_size = 0
    if compressed
        max_zoom_uncomp_size = items_per_slot * sizeof(BigBedZoomData)
        uncompress_buf_size = max(max_block_size, max_zoom_uncomp_size)
    end

    # go back and rewrite header
    seek(out, 0)
    write(out,
        BigBedHeader(fmt == BigBed ? BIGBED_MAGIC : BIGWIG_MAGIC,
                     BIGBED_CURRENT_VERSION, zoom_levels,
                     chrom_tree_offset, data_offset, index_offset,
                     field_count, min(field_count, 12), as_offset,
                     total_summary_offset, uncompress_buf_size, 0))

    # write summary headers with data
    for i in 1:zoom_levels
        write(out, BigBedZoomHeader(zoom_amounts[i], 0, zoom_data_offsets[i],
                                    zoom_index_offsets[i]))
    end

    # write rest of summary header with no data
    for i in zoom_levels:(BIGBED_MAX_ZOOM_LEVELS-1)
        write(out, BigBedZoomHeader(0, 0, 0, 0))
    end

    # write total summary
    seek(out, total_summary_offset)
    write(out, total_summary)

    # write end signature
    seekend(out)
    write(out, convert(UInt32, fmt == BigBed ? BIGBED_MAGIC : BIGWIG_MAGIC))
end
