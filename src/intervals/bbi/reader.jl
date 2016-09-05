# Reader
# ======

"""
    BigBedReader(input::IO)

Create a data reader of the BigBed file format.

# Arguments
* `input`: data source
"""
type BigBedReader <: Bio.IO.AbstractReader
    stream::BufferedInputStream
    header::BigBedHeader
    zoom_headers::Vector{BigBedZoomHeader}
    autosql::String
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

function Base.iteratorsize(::BigBedReader)
    return Base.SizeUnknown()
end

function Base.eltype(::Type{BigBedReader})
    return BEDInterval
end

function Bio.IO.stream(reader::BigBedReader)
    return reader.stream
end

function BigBedReader(input::IO)
    return init_bigbed_reader(input)
end

# Initialize BgiBedData for reading.  Once opened, entries can be read from the
# file either by iterating over it, or by indexing into it with an interval.
function init_bigbed_reader(stream::IO)
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

    return BigBedReader(BufferedInputStream(stream), header, zoom_headers, autosql, summary,
                      btree_header, rtree_header, data_count,
                      BigBedBTreeInternalNode[BigBedBTreeInternalNode(btree_header.key_size)
                                              for _ in 1:btree_header.block_size],
                      BigBedBTreeLeafNode[BigBedBTreeLeafNode(btree_header.key_size)
                                          for _ in 1:btree_header.block_size],
                      Array(UInt8, btree_header.key_size),
                      Array(Vector{UInt8}, btree_header.block_size),
                      Array(UInt8, header.uncompress_buf_size))
end

# Return all sequence (name, id, size) tuples in a BigBed B-tree.
function first_btree_leaf_position(bb::BigBedReader)
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


# Parser
# ------

type BigBedDataReader <: Bio.IO.AbstractReader
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

    function BigBedDataReader(input::BufferedInputStream;
                              seq_names::Nullable{Vector{StringField}}=Nullable{Vector{StringField}}(),
                              assumed_seqname::Nullable{StringField}=Nullable{StringField}())
        cs = _bigbedparser_start
        return new(Ragel.State(cs, input), 0, 0.0, 0.0, 0.0, 1, 1, seq_names, assumed_seqname)
    end
end

function Bio.IO.stream(reader::BigBedDataReader)
    return reader.state.stream
end

include("parser.jl")


# Iterator
# --------

# An iterator over all entries in a BigBed file.
type BigBedIteratorState
    seq_names::Vector{StringField}
    data_count::Int
    data_num::Int
    data_offset::UInt
    reader::BigBedDataReader
    reader_isdone::Bool
    next_interval::Interval{BEDMetadata}
end

function Base.start(bb::BigBedReader)
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

    reader = BigBedDataReader(
        BufferedInputStream(bb.uncompressed_data, unc_block_size),
        seq_names=Nullable(seq_names))

    next_interval = BEDInterval()
    reader_isdone = isnull(tryread!(reader, next_interval))

    return BigBedIteratorState(seq_names, data_count, 1,
                               zlib_stream.source.zstream.total_in,
                               reader, reader_isdone, next_interval)
end

function Base.next(bb::BigBedReader, state::BigBedIteratorState)
    value = copy(state.next_interval)

    state.data_num += 1
    if state.data_num > state.data_count
        return value, state
    end

    state.reader_isdone = isnull(tryread!(state.reader, state.next_interval))

    if state.reader_isdone
        seek(bb.stream, bb.header.full_data_offset + state.data_offset + sizeof(UInt64))
        zlib_stream = ZlibInflateInputStream(bb.stream, reset_on_end=false)

        unc_block_size = readbytes!(zlib_stream, bb.uncompressed_data,
                                    length(bb.uncompressed_data))
        state.reader = BigBedDataReader(
             BufferedInputStream(bb.uncompressed_data, unc_block_size),
             seq_names=Nullable(state.seq_names))
        state.data_offset += zlib_stream.source.zstream.total_in

        state.reader_isdone = isnull(tryread!(state.reader, state.next_interval))
        @assert !state.reader_isdone
    end

    return value, state
end

function Base.done(bb::BigBedReader, state::BigBedIteratorState)
    return state.data_num > state.data_count
end
