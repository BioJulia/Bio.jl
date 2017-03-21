# Intersection
# ============

# An iterator over entries in a BigBed file that intersect a given interval.
# Constructed by indexing into `BigBedReader` with an interval.
type BigBedIntersectIterator
    bb::BigBedReader

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

    reader::Nullable{BigBedDataReader}
    nextinterval::BEDInterval
    done::Bool
end

function Base.iteratorsize(::Type{BigBedIntersectIterator})
    return Base.SizeUnknown()
end

function eachoverlap(bb::BigBedReader, query::Interval)
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
                                   Nullable{BigBedDataReader}(),
                                   BEDInterval(), false)
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

# Find the given seqname in the BigBed file's index and read the corresponding
# sequence id and length.
function lookup_seqname(bb::BigBedReader, seqname::AbstractString)
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

function find_next_intersection!(it::BigBedIntersectIterator)
    it.done = true
    while !isnull(it.reader) || it.block_num < length(it.blocks)
        it.done = true
        while !isnull(it.reader)
            reader = get(it.reader)
            # Run the reader until we find an intersection
            it.done = isnull(tryread!(reader, it.nextinterval))
            if it.done
                it.reader = Nullable{BigBedDataReader}()
                break
            end

            if reader.chrom_id == it.query_chrom_id &&
               it.nextinterval.first <= it.query_last &&
               it.nextinterval.last >= it.query_first
                it.done = false
                return
            end
        end

        if it.block_num < length(it.blocks)
            # advance to the next block of interest ant initialize a new reader
            it.block_num += 1
            block_offset, block_size = it.blocks[it.block_num]

            seek(it.bb.stream, block_offset)
            @assert block_size <= length(it.bb.uncompressed_data)
            unc_block_size = readbytes!(Libz.ZlibInflateInputStream(it.bb.stream, reset_on_end=false),
                                        it.bb.uncompressed_data,
                                        length(it.bb.uncompressed_data))
            it.reader = BigBedDataReader(
                BufferedInputStream(it.bb.uncompressed_data, unc_block_size),
                assumed_seqname=Nullable(it.query_seqname))
        else
            it.done = true
            return
        end
    end
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
