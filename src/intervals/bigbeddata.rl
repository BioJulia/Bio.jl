
type BigBedData <: IntervalStream{BEDMetadata}
    reader::BufferedReader
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


module BigBedDataParserImpl

import Bio.Ragel, Zlib
import Bio.Intervals: Strand, STRAND_NA, BEDInterval, BEDMetadata

using Color, Switch

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
        input.chrom_id = unsafe_load(convert(Ptr{UInt32}, pointer(state.reader.buffer, m)))
    }

    action first {
        m = Ragel.@unmark!
        input.first = unsafe_load(convert(Ptr{UInt32}, pointer(state.reader.buffer, m)))
    }

    action last {
        m = Ragel.@unmark!
        input.last = unsafe_load(convert(Ptr{UInt32}, pointer(state.reader.buffer, m)))
    }

    action name        { input.name         = Nullable{AbstractString}(Ragel.@asciistring_from_mark!) }
    action score       { input.score        = Ragel.@int64_from_mark! }
    action strand      { input.strand       = convert(Strand, Ragel.@char) }
    action thick_first { input.thick_first  = (Ragel.@int64_from_mark!) + 1 }
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
        push!(get(input.block_firsts), (Ragel.@int64_from_mark!) + 1)
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

    chrom_id::UInt32
    first::Int64
    last::Int64
    strand::Strand

    red::Float32
    green::Float32
    blue::Float32
    name::Nullable{AbstractString}
    score::Nullable{Int}
    thick_first::Nullable{Int}
    thick_last::Nullable{Int};
    item_rgb::Nullable{RGB{Float32}};
    block_count::Nullable{Int}
    block_sizes::Nullable{Vector{Int}}
    block_firsts::Nullable{Vector{Int}}

    function BigBedDataParser(input::Vector{UInt8}, len::Integer)
        %%write init;

        return new(Ragel.State(cs, input, false, len),
                   0, 0, 0, STRAND_NA, 0.0, 0.0, 0.0,
                   Nullable{AbstractString}(), Nullable{Int}(), Nullable{Int}(),
                   Nullable{Int}(), Nullable{RGB{Float32}}(), Nullable{Int}(),
                   Nullable{Vector{Int}}(), Nullable{Vector{Int}}())
    end
end


function takevalue!(input::BigBedDataParser, seqname::ASCIIString)

    interval = BEDInterval(
        seqname,
        input.first + 1, input.last, input.strand,
        BEDMetadata(input.name, input.score, input.thick_first,
                    input.thick_last, input.item_rgb,
                    input.block_count, input.block_sizes,
                    input.block_firsts))

    input.strand = STRAND_NA
    input.name = Nullable{ASCIIString}()
    input.score = Nullable{Int}()
    input.thick_first = Nullable{Int}()
    input.thick_last = Nullable{Int}()
    input.item_rgb = Nullable{RGB{Float32}}()
    input.block_count = Nullable{Int}()
    input.block_sizes = Nullable{Vector{Int}}()
    input.block_firsts = Nullable{Vector{Int}}()

    return interval
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
