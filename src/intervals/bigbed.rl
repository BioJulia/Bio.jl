

# The BigBed format is documented in
#   Kent, W. James, et al. "BigWig and BigBed: # enabling browsing of large
#   distributed datasets." Bioinformatics 26.17 (2010): 2204-2207.
#
# The low level details are documented in a series of tables in the supplement
# of that paper. The immutable types defined below are exactly the data layout
# described in those tables. They are labeled with the table they correspand to.


const BIGBED_MAGIC = 0x8789F2EB
const BIGWIG_MAGIC = 0x888FFC26
const BIGBED_BTREE_MAGIC = 0x78CA8C91
const BIGBED_RTREE_MAGIC = 0x2468ACE0


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


# Supplemental Table 14
immutable BigBedRTreeHeader
    magic::Uint32
    block_size::Uint32
    item_count::Uint64
    start_chrom_ix::Uint32
    start_base::Uint32
    end_chrom_ix::Uint32
    end_base::Uint32
    end_file_offset::Uint64
    items_per_slot::Uint32
    reserved::Uint32
end


function read(io::IO, ::Type{BigBedRTreeHeader})
    return BigBedRTreeHeader(
        read(io, Uint32), read(io, Uint32), read(io, Uint64), read(io, Uint32),
        read(io, Uint32), read(io, Uint32), read(io, Uint32), read(io, Uint64),
        read(io, Uint32), read(io, Uint32))
end


# Supplemental Table 15
immutable BigBedRTreeNode
    isleaf::Uint8
    reserved::Uint8
    count::Uint16
end


function read(io::IO, ::Type{BigBedRTreeNode})
    return BigBedRTreeNode(read(io, Uint8), read(io, Uint8), read(io, Uint16))
end


# Supplemental Table 16
immutable BigBedRTreeLeafNode
    start_chrom_ix::Uint32
    start_base::Uint32
    end_chrom_ix::Uint32
    end_base::Uint32
    data_offset::Uint64
    data_size::Uint64
end


function read(io::IO, ::Type{BigBedRTreeLeafNode})
    return BigBedRTreeLeafNode(
        read(io, Uint32), read(io, Uint32), read(io, Uint32), read(io, Uint32),
        read(io, Uint64), read(io, Uint64))
end


# Supplemental Table 17
immutable BigBedRTreeInternalNode
    start_chrom_ix::Uint32
    start_base::Uint32
    end_chrom_ix::Uint32
    end_base::Uint32
    data_offset::Uint64
end


function read(io::IO, ::Type{BigBedRTreeInternalNode})
    return BigBedRTreeInternalNode(
        read(io, Uint32), read(io, Uint32), read(io, Uint32), read(io, Uint32),
        read(io, Uint64))
end


using Bio: BufferedReader

immutable BigBed <: FileFormat end


type BigBedData
    reader::BufferedReader
    header::BigBedHeader
    zoom_headers::Vector{BigBedZoomHeader}
    autosql::String
    summary::BigBedTotalSummary
    btree_header::BigBedBTreeHeader
    rtree_header::BigBedRTreeHeader
    data_count::Uint32

    # preallocated space for reading and searchig the B-tree
    btree_internal_nodes::Vector{BigBedBTreeInternalNode}
    btree_leaf_nodes::Vector{BigBedBTreeLeafNode}
    key::Vector{Uint8}
    node_keys::Vector{Vector{Uint8}}
    uncompressed_data::Vector{Uint8}
end


module BigBedDataParserImpl

import Bio.Ragel, Zlib
import Bio.Intervals: Strand, STRAND_NA, BEDInterval, BEDMetadata
using Color, Compat, Switch

# Parser for data blocks in a BigBed file. This is very similar
# to the BED parser in bed.rl, with the following exceptions:
#
#    * BigBed has binary chrom_index, start, and end, insteado of ASCII
#      chromosome name, start, end.
#    * BigBed entries are null ('\0') terminated, rather than newline separated.
#
%%{
    machine bigbed;

    action yield {
        yield = true
        # // fbreak causes will cause the pushmark action for the next seqname
        # // to be skipped, so we do it here
        Ragel.@mark!
        fbreak;
    }

    action mark { Ragel.@mark! }

    action chrom_id {
        m = Ragel.@unmark!
        input.chrom_id = unsafe_load(convert(Ptr{Uint32}, pointer(state.reader.buffer, m)))
    }

    action first {
        m = Ragel.@unmark!
        input.first = unsafe_load(convert(Ptr{Uint32}, pointer(state.reader.buffer, m))) - 1
    }

    action last {
        m = Ragel.@unmark!
        input.last = unsafe_load(convert(Ptr{Uint32}, pointer(state.reader.buffer, m)))
    }

    action name        { input.name         = Nullable{String}(Ragel.@asciistring_from_mark!) }
    action score       { input.score        = Ragel.@int64_from_mark! }
    action strand      { input.strand       = convert(Strand, Ragel.@char) }
    action thick_first { input.thick_first  = Ragel.@int64_from_mark! }
    action thick_last  { input.thick_last   = Ragel.@int64_from_mark! }
    action item_rgb_r  { input.red = input.green = input.blue = (Ragel.@int64_from_mark!) / 255.0 }
    action item_rgb_g  { input.green        = (Ragel.@int64_from_mark!) / 255.0 }
    action item_rgb_b  { input.blue         = (Ragel.@int64_from_mark!) / 255.0 }
    action item_rgb    { input.item_rgb     = RGB{Float32}(input.red, input.green, input.blue ) }
    action block_count { input.block_count  = Ragel.@int64_from_mark! }
    action block_size {
        if isnull(input.block_sizes)
            input.block_sizes = Array(Int, 0)
        end
        push!(get(input.block_sizes), (Ragel.@int64_from_mark!))
    }
    action block_first {
        if isnull(input.block_firsts)
            input.block_firsts = Array(Int, 0)
        end
        push!(get(input.block_firsts), (Ragel.@int64_from_mark!))
    }

    hspace       = [ \t\v];

    chrom_id     = any{4}   >mark     %chrom_id;
    first        = any{4}   >mark     %first;
    last         = any{4}   >mark     %last;
    name         = [ -~]*   >mark     %name;
    score        = digit+   >mark     %score;
    strand       = [+\-\.?] >strand;
    thick_first  = digit+   >mark     %thick_first;
    thick_last   = digit+   >mark     %thick_last;

    item_rgb_r   = digit+   >mark     %item_rgb_r;
    item_rgb_g   = digit+   >mark     %item_rgb_g;
    item_rgb_b   = digit+   >mark     %item_rgb_b;
    item_rgb     = item_rgb_r (hspace* ',' hspace* item_rgb_g hspace* ',' hspace* item_rgb_b)? %item_rgb;

    block_count  = digit+   >mark    %block_count;

    block_size   = digit+   >mark    %block_size;
    block_sizes  = block_size (',' block_size)* ','?;

    block_first  = digit+   >mark    %block_first;
    block_firsts = block_first (',' block_first)* ','?;

    bed_entry = chrom_id first last (
                    name ( '\t' score ( '\t' strand ( '\t' thick_first (
                    '\t' thick_last ( '\t' item_rgb ( '\t' block_count (
                    '\t' block_sizes ( '\t' block_firsts )? )? )? )? )? )? )? )?
                    )? '\0';
    main := (bed_entry %yield)*;
}%%


%%write data;

type BigBedDataParser
    state::Ragel.State

    chrom_id::Uint32
    first::Int64
    last::Int64
    strand::Strand

    red::Float32
    green::Float32
    blue::Float32
    name::Nullable{String}
    score::Nullable{Int}
    thick_first::Nullable{Int}
    thick_last::Nullable{Int};
    item_rgb::Nullable{RGB{Float32}};
    block_count::Nullable{Int}
    block_sizes::Nullable{Vector{Int}}
    block_firsts::Nullable{Vector{Int}}

    function BigBedDataParser(input::Vector{Uint8}, len::Integer)
        %%write init;

        return new(Ragel.State(cs, input, len),
                   0, 0, 0, STRAND_NA, 0.0, 0.0, 0.0,
                   Nullable{String}(), Nullable{Int}(), Nullable{Int}(),
                   Nullable{Int}(), Nullable{RGB{Float32}}(), Nullable{Int}(),
                   Nullable{Vector{Int}}(), Nullable{Vector{Int}}())
    end
end



Ragel.@generate_read_fuction("bigbed", BigBedDataParser, BigBedDataEntry,
    begin
        @inbounds begin
            %%write exec;
        end
    end,
    begin
    end)


end # module BigBedDataParserImpl



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
    btree_header = read(reader, BigBedBTreeHeader)
    if btree_header.magic != BIGBED_BTREE_MAGIC
        error("BigBed B-Tree magic number was incorrect. File may be corrupt or malformed.")
    end

    # data_count
    seek(reader, header.full_data_offset + 1)
    data_count = read(reader, Uint32)

    # r-tree header
    seek(reader, header.full_index_offset + 1)
    rtree_header = read(reader, BigBedRTreeHeader)
    if rtree_header.magic != BIGBED_RTREE_MAGIC
        error("BigBed R-Tree magic number was incorrect. File may be corrupt or malformed.")
    end

    return BigBedData(reader, header, zoom_headers, autosql, summary,
                      btree_header, rtree_header, data_count,
                      BigBedBTreeInternalNode[BigBedBTreeInternalNode(btree_header.key_size)
                                              for _ in 1:btree_header.block_size],
                      BigBedBTreeLeafNode[BigBedBTreeLeafNode(btree_header.key_size)
                                          for _ in 1:btree_header.block_size],
                      Array(Uint8, btree_header.key_size),
                      Array(Vector{Uint8}, btree_header.block_size),
                      Array(Uint8, header.uncompressed_buf_size))
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
An iterator over entries in a BigBed file that intersect a given interval.

Constructed by indexing into `BigBedData` with an interval.
""" ->
type BigBedIntersectIterator
    bb::BigBedData

    query_seqname::String
    query_first::Int64
    query_last::Int64
    query_chrom_id::Uint32
    query_chrom_size::Uint32

    # (offset, size) pairs giving the extents of data blocks that may contain
    # (overlapping intervals
    blocks::Vector{(@compat Tuple{Uint64, Uint64})}

    # Index of block currently being parsed
    block_num::Int

    parser::Nullable{BigBedDataParserImpl.BigBedDataParser}
    nextinterval::Nullable{BEDInterval}
end


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
            for i in 1:bb.btree_header.block_size
                read!(bb.reader, bb.btree_internal_nodes[i])
                bb.node_keys[i] = bb.btree_internal_nodes[i].key
            end
            i = searchsortedfirst(bb.node_keys, bb.key, Base.Order.Lt(memisless))
            if !(1 <= i <= bb.btree_header.block_size)
                break
            end

            seek(bb.reader, bb.btree_internal_nodes[i].child_offset + 1)
        else
            for i in 1:bb.btree_header.block_size
                read!(bb.reader, bb.btree_leaf_nodes[i])
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


function getindex(bb::BigBedData, query::Interval)
    seqname, first, last = query.seqname, query.first, query.last
    chrom_id, chrom_size = lookup_seqname(bb, seqname)

    # stack of nodes yet to visit
    root_position = bb.header.full_index_offset + sizeof(BigBedRTreeHeader) + 1
    node_positions = Uint64[root_position]

    # stack of data blocks that may contain overlapping intervals
    blocks = (@compat Tuple{Uint64, Uint64})[]

    seek(bb.reader, bb.header.full_index_offset + sizeof(BigBedRTreeHeader) + 1)
    while !isempty(node_positions)
        pos = pop!(node_positions)
        seek(bb.reader, pos)
        node = read(bb.reader, BigBedRTreeNode)
        if node.isleaf == 0
            for i in 1:node.count
                internal_node = read(bb.reader, BigBedRTreeInternalNode)
                if internal_node.start_chrom_ix <= chrom_id <= internal_node.end_chrom_ix &&
                    (chrom_id < internal_node.end_chrom_ix || first <= internal_node.end_base) &&
                    (chrom_id > internal_node.start_chrom_ix || last - 1 >= internal_node.start_base)
                    push!(node_positions, internal_node.data_offset + 1)
                end
            end
        else
            for i in 1:node.count
                leaf_node = read(bb.reader, BigBedRTreeLeafNode)
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
                                   Nullable{BigBedDataParserImpl.BigBedDataParser}(),
                                   Nullable{BEDInterval}())
end


function find_next_intersection!(it::BigBedIntersectIterator)
    it.nextinterval = Nullable{BEDInterval}()
    while !isnull(it.parser) || it.block_num < length(it.blocks)
        while !isnull(it.parser)
            parser = get(it.parser)
            # Run the parser until we find an intersection
            isdone = !BigBedDataParserImpl.advance!(parser)
            if isdone
                it.parser = Nullable{BigBedDataParserImpl.BigBedDataParser}()
                break
            end

            # Note: (parser.first, parser.last) is 0-based, end-exclusive,
            # while (it.query_first, it.query_last) is 1-based, end-inclusive.
            if parser.chrom_id == it.query_chrom_id &&
               parser.first <= it.query_last &&
               parser.last - 1 >= it.query_first
                it.nextinterval = BEDInterval(
                    it.query_seqname, parser.first + 1, parser.last,
                    parser.strand,
                    BEDMetadata(parser.name, parser.score, parser.thick_first,
                                parser.thick_last, parser.item_rgb,
                                parser.block_count, parser.block_sizes,
                                parser.block_firsts))
                return
            end
        end

        if it.block_num < length(it.blocks)
            # advance to the next block of interest ant initialize a new parser
            it.block_num += 1
            block_offset, block_size = it.blocks[it.block_num]

            seek(it.bb.reader, block_offset + 1)
            @assert block_size <= length(it.bb.uncompressed_data)
            unc_block_size = readbytes!(Zlib.Reader(it.bb.reader), it.bb.uncompressed_data,
                                        length(it.bb.uncompressed_data))
            it.parser = BigBedDataParserImpl.BigBedDataParser(
                it.bb.uncompressed_data, unc_block_size)
        else
            break
        end
    end
end


function start(it::BigBedIntersectIterator)
    find_next_intersection!(it)
    nothing
end


function next(it::BigBedIntersectIterator, ::Nothing)
    value = it.nextinterval
    find_next_intersection!(it)
    return value, nothing
end


function done(it::BigBedIntersectIterator, ::Nothing)
    return isnull(it.nextinterval)
end


