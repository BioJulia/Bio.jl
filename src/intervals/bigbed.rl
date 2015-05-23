

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

const BIGBED_MAX_ZOOM_LEVELS = 10
const BIGBED_CURRENT_VERSION = 4


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
    uncompress_buf_size::Uint32
    reserved::Uint64
end


function read(io::IO, ::Type{BigBedHeader})
    return BigBedHeader(
        read(io, Uint32), read(io, Uint16), read(io, Uint16),
        read(io, Uint64), read(io, Uint64), read(io, Uint64),
        read(io, Uint16), read(io, Uint16), read(io, Uint64),
        read(io, Uint64), read(io, Uint32), read(io, Uint64))
end


function write(io::IO, header::BigBedHeader)
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
    reduction_level::Uint32
    reserved::Uint32
    data_offset::Uint64
    index_offset::Uint64
end


function read(io::IO, ::Type{BigBedZoomHeader})
    return BigBedZoomHeader(
        read(io, Uint32), read(io, Uint32), read(io, Uint64), read(io, Uint64))
end


function write(io::IO, header::BigBedZoomHeader)
    write(io, header.reduction_level)
    write(io, header.reserved)
    write(io, header.data_offset)
    write(io, header.index_offset)
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


# Supplemental Table 19
immutable BigBedZoomData
    chrom_id::Uint32
    chrom_start::Uint32
    chrom_end::Uint32
    valid_count::Uint32
    min_val::Float32
    max_val::Float32
    sum_data::Float32
    sum_squares::Float32
end


function read(io::IO, ::Type{BigBedZoomData})
    return BigBedZoomData(
        read(io, Uint32), read(io, Uint32), read(io, Uint32), read(io, Uint32),
        read(io, Float32), read(io, Float32), read(Float32), read(Float32))
end


function write(io::IO, data::BigBedZoomData)
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


function write(io::IO, header::BigBedBTreeHeader)
    write(io, header.magic, header.block_size, header.key_size,
          header.val_size, header.item_count, header.reserved)
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


function write(io::IO, node::BigBedBTreeNode)
    write(io, node.isleaf)
    write(io, node.reserved)
    write(io, node.count)
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


function write(io::IO, header::BigBedRTreeHeader)
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
    isleaf::Uint8
    reserved::Uint8
    count::Uint16
end


function read(io::IO, ::Type{BigBedRTreeNode})
    return BigBedRTreeNode(read(io, Uint8), read(io, Uint8), read(io, Uint16))
end


function write(io::IO, node::BigBedRTreeNode)
    write(io, node.isleaf)
    write(io, node.reserved)
    write(io, node.count)
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


function write(io::IO, node::BigBedRTreeLeafNode)
    write(io, node.start_chrom_ix)
    write(io, node.start_base)
    write(io, node.end_chrom_ix)
    write(io, node.end_base)
    write(io, node.data_offset)
    write(io, node.data_size)
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


function write(io::IO, node::BigBedRTreeInternalNode)
    write(io, node.start_chrom_ix)
    write(io, node.start_base)
    write(io, node.end_chrom_ix)
    write(io, node.end_base)
    write(io, node.data_offset)
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

    # preallocated space for reading and searching the B-tree
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
    if header.auto_sql_offset == 0
        autosql = ""
    else
        seek(reader, header.auto_sql_offset + 1)
        autosql_buf = IOBuffer()
        while (c = read(reader, Uint8)) != 0x00
            write(autosql_buf, c)
        end
        # TODO: eventually we should parse this and do something useful with it
        autosql = takebuf_string(autosql_buf)
    end

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
                      Array(Uint8, header.uncompress_buf_size))
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


# BigBed Output
# -------------


function write_zeros(out::IO, n::Integer)
    for i in 1:n
        write(out, 0x00)
    end
end


immutable BigBedChromInfo
    name::ASCIIString   # chromosome name
    id::Uint32          # unique id for chromosome
    size::Uint32        # size of chromosome
    item_count::Uint64  # number of items
end


immutable BigBedBounds
    offset::Uint64
    chrom_ix::Uint32
    start::Uint32
    stop::Uint32
end


function bigbed_chrom_info(intervals::IntervalCollection, chrom_sizes::Dict)
    info = Array(BigBedChromInfo, length(intervals.trees))
    seqnames = collect(ASCIIString, keys(intervals.trees))
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
    bytes_in_index_block = block_header_size + block_size * (key_size * 8)
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
        endix = min(length(items), i + node_size_per)
        for j in i:slot_size_per:end_ix
            item = items[j]
            write(out, item.name)
            # pad name
            for k in length(item.name)+1:key_size
                write(out, 0x00)
            end
            write(out, convert(Uint64, next_child))
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
            write(out, convert(Uint32, item.id))
            write(out, convert(Uint32, item.size))
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
                                 chrom_sizes::Dict=Dict{String, Int64}())
    # See bbiWriteChromInfo in bbiWrite.c
    length(intervals.trees)
    for seqname in keys(intervals.trees)
        if !in(seqname, chrom_sizes)
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

                write(buf, convert(Uint32, info.id))
                write(buf, convert(Uint32, interval_start))
                write(buf, convert(Uint32, interval_end))
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
                writer = Zlib.Writer(out, 6)
                write(writer, block)
                close(writer)
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


type BigBedRTree
    next::Nullable{BigBedRTree}
    children::Nullable{BigBedRTree}
    parent::Nullable{BigBedRTree}
    start_chrom_ix::Uint32
    start_base::Uint32
    end_chrom_ix::Uint32
    end_base::Uint32
    start_file_offset::Uint64
    end_file_offset::Uint64

    function BigBedRTree(next, children, parent, start_chrom_ix, start_base,
                         end_chrom_ix, end_base, start_file_offset, end_file_offset)
        return new(next, children, parent, start_chrom_ix, end_chrom_ix,
                   start_base, end_base, start_file_offset, end_file_offset)
    end

    function BigBedRTree()
        return new()
    end
end


function copy(tree::BigBedRTree)
    return BigBedRTree(tree.next, tree.children, tree.parent,
                       tree.start_chrom_ix, tree.start_base,
                       tree.end_chrom_ix, tree.end_base,
                       tree.start_file_offset, tree.end_file_offset)
end


function length(list::Nullable{BigBedRTree})
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


function write(out::IO, tree::BigBedRTree, block_size, level_count)
    # See writeTreeToOpenFile in cirTree.c

    level_sizes = zeros(Int, level_count)
    bigbed_calc_rtree_level_sizes!(level_sizes, Nullable{BigBedRTree}(tree), 1, level_count)

    level_offsets = Array(Uint64, level_count)

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


function reverse!(tree::Nullable{BigBedRTree})
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


# Write out data reduced by factor of initialReduction.  Also calculate and keep
# in memory next reduction level.  This is more work than some ways, but it
# keeps us from having to keep the first reduction entirely in memory.
function bigbed_write_reduced_once_return_reduced_twice(
                out::IO, chrom_info::Vector{BigBedChromInfo},
                field_count, intervals::IntervalCollection,
                initial_reduction, initial_reduced_count,
                zoom_increment, block_size, items_per_slot, compressed)

    double_reduction_size = initial_reduction * zoom_increment
    bounds_pt = Array(BigBedBounds, initial_reduced_count)

    ret_data_start = position(out)
    write(out, convert(Uint32, initial_reduced_count))

    # This gets a little complicated I'm afraid.  The strategy is to:
    #  1) Build up a range tree that represents coverage depth on that chromosome
    #     This also has the nice side effect of getting rid of overlaps.
    #  2) Stream through the range tree, outputting the initial summary level and
    #     further reducing. 
    #
    # TODO
end


function bed_field_count(intervals::IntervalCollection{BEDMetadata})
    num_fields = 0
    for interval in intervals
        interval_num_fields = 3

        if !isnull(interval.metadata.name)
            interval_num_fields += 1
        else @goto finished end

        if !isnull(interval.metadata.score)
            interval_num_fields += 1
        else @goto finished end

        if interval.strand != STRAND_NA
            interval_num_fields += 1
        else @goto finished end

        if !isnull(interval.metadata.thick_first)
            interval_num_fields += 1
        else @goto finished end

        if !isnull(interval.metadata.thick_last)
            interval_num_fields += 1
        else @goto finished end

        if !isnull(interval.metadata.item_rgb)
            interval_num_fields += 1
        else @goto finished end

        if !isnull(interval.metadata.block_count)
            interval_num_fields += 1
        else @goto finished end

        if !isnull(interval.metadata.block_sizes)
            interval_num_fields += 1
        else @goto finished end

        if !isnull(interval.metadata.block_firsts)
            interval_num_fields += 1
        else @goto finished end

        @label finished

        if num_fields == 0 || num_fields > interval_num_fields
            num_fields = interval_num_fields
        end
    end
    return num_fields
end


function write(out::IO, ::Type{BigBed}, intervals::IntervalCollection;
               block_size::Int=256, items_per_slot::Int=512,
               compressed::Bool=true)
    # See the function bbFileCreate in bedToBigBed.c in kent to see the only
    # other implementation of this I'm aware of.

    if !applicable(seek, out, 1)
        error("BigBed can only be written to seekable output streams.")
    end

    field_count = bed_field_count(intervals)

    # write dummy headers, come back and fill them later
    write_zeros(out, sizeof(BigBedHeader))
    write_zeros(out, BIGBED_MAX_ZOOM_LEVELS * sizeof(BigBedZoomHeader))

    # TODO: optional autoSql specification
    as_offset = 0

    total_summary_offset = position(out)
    write_zeros(out, sizeof(BigBedTotalSummary))

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
    res = max(ave_span, min_zoom)
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
    write(out, convert(Uint64, length(intervals)))
    block_count = bigbed_count_sections_needed(chrom_info, items_per_slot)
    bounds = Array(BigBedBounds, block_count)
    max_block_size = bigbed_write_blocks(out, intervals, chrom_info, items_per_slot,
                                         bounds, block_count, compressed,
                                         res_try_count, res_scales, res_sizes)

    # write out primary data index
    index_offset = position(out)
    bigbed_write_index(out, bounds, block_size, 1, index_offset)

    # declare arrays and vars that track the zoom levels we actually output
    # TODO: build and write zoom level data
    zoom_levels = 0
    if false
        zoom_amounts = Array(Uint32, BIGBED_MAX_ZOOM_LEVELS)
        zoom_data_offsets = Array(Uint64, BIGBED_MAX_ZOOM_LEVELS)
        zoom_index_offsets = Array(Uint64, BIGBED_MAX_ZOOM_LEVELS)
        zoom_levels = 0

        @show ave_span
        @show res_try_count
        @show res_sizes

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

            @show initial_reduction
            @show initial_reduced_count

            if initial_reduction > 0
                zoom_increment = 4
                reoomed_list = bigbed_write_reduced_once_return_reduced_twice()
                zoom_amounts[1] = initial_reduction
                zoom_levels = 1

                while zoom_levels < BIGBED_MAX_ZOOM_LEVELS
                    # TODO
                end
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
        BigBedHeader(BIGBED_MAGIC, BIGBED_CURRENT_VERSION, zoom_levels,
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
    # TODO: summary

    # write end signature
    seekend(out)
    write(out, convert(Uint32, BIGBED_MAGIC))
end


