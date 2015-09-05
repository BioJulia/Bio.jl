
type BigBedData <: IntervalStream{BEDMetadata}
    stream::BufferedInputStream
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

import Bio.Ragel
using Bio: StringField
using Bio.Intervals: Strand, STRAND_NA, BEDInterval, BEDMetadata
using Colors, Compat, Switch, BufferedStreams

# Parser for data blocks in a BigBed file. This is very similar
# to the BED parser in bed.rl, with the following exceptions:
#
#    * BigBed has binary chrom_index, start, and end, insteado of ASCII
#      chromosome name, start, end.
#    * BigBed entries are null ('\0') terminated, rather than newline separated.
#
%%{
    machine bigbed;

    action finish_match {
        input.block_size_idx = 1
        input.block_first_idx = 1

        yield = true
        # // fbreak causes will cause the pushmark action for the next seqname
        # // to be skipped, so we do it here
        Ragel.@anchor!
        fbreak;
    }

    action anchor { Ragel.@anchor! }
    action optional_field { output.metadata.used_fields += 1 }

    action chrom_id {
        input.chrom_id = Ragel.@load_from_anchor!(UInt32)
        if !isnull(input.seq_names)
            output.seqname = copy(get(input.seq_names)[input.chrom_id + 1])
        else
            output.seqname = copy(get(input.assumed_seqname))
        end
        output.metadata.used_fields = 0
    }

    action first       { output.first = 1 + Ragel.@load_from_anchor!(UInt32) }
    action last        { output.last = Ragel.@load_from_anchor!(UInt32) }
    action name        { Ragel.@copy_from_anchor!(output.metadata.name) }
    action score       { output.metadata.score = Ragel.@int64_from_anchor! }
    action strand      { output.strand = convert(Strand, (Ragel.@char)) }
    action thick_first { output.metadata.thick_first = 1 + Ragel.@int64_from_anchor! }
    action thick_last  { output.metadata.thick_last = Ragel.@int64_from_anchor!  }
    action item_rgb_r  { input.red = input.green = input.blue = (Ragel.@int64_from_anchor!) / 255.0 }
    action item_rgb_g  { input.green = (Ragel.@int64_from_anchor!) / 255.0 }
    action item_rgb_b  { input.blue = (Ragel.@int64_from_anchor!) / 255.0 }
    action item_rgb    { output.metadata.item_rgb = RGB{Float32}(input.red, input.green, input.blue) }
    action block_count {
        output.metadata.block_count = Ragel.@int64_from_anchor!

        if (output.metadata.block_count > length(output.metadata.block_sizes))
            resize!(output.metadata.block_sizes, output.metadata.block_count)
        end

        if (output.metadata.block_count > length(output.metadata.block_firsts))
            resize!(output.metadata.block_firsts, output.metadata.block_count)
        end
    }

    action block_size {
        if input.block_size_idx > length(output.metadata.block_sizes)
            error("More size blocks encountered than BED block count field suggested.")
        end
        output.metadata.block_sizes[input.block_size_idx] = Ragel.@int64_from_anchor!
        input.block_size_idx += 1
    }

    action block_first {
        if input.block_first_idx > length(output.metadata.block_firsts)
            error("More start blocks encountered than BED block count field suggested.")
        end
        output.metadata.block_firsts[input.block_first_idx] = 1 + Ragel.@int64_from_anchor!
        input.block_first_idx += 1
    }

    hspace       = [ \t\v];

    chrom_id     = any{4}   >anchor     %chrom_id;
    first        = any{4}   >anchor     %first;
    last         = any{4}   >anchor     %last;
    name         = [ -~]+   >anchor     %name  %optional_field;
    score        = digit+   >anchor     %score %optional_field;
    strand       = [+\-\.?] >strand %optional_field;
    thick_first  = digit+   >anchor     %thick_first %optional_field;
    thick_last   = digit+   >anchor     %thick_last  %optional_field;

    item_rgb_r   = digit+   >anchor     %item_rgb_r;
    item_rgb_g   = digit+   >anchor     %item_rgb_g;
    item_rgb_b   = digit+   >anchor     %item_rgb_b;
    item_rgb     = item_rgb_r (hspace* ',' hspace* item_rgb_g hspace* ',' hspace* item_rgb_b)? %item_rgb %optional_field;

    block_count  = digit+   >anchor    %block_count %optional_field;

    block_size   = digit+   >anchor    %block_size;
    block_sizes  = block_size (',' block_size)* ','? %optional_field;

    block_first  = digit+   >anchor    %block_first;
    block_firsts = block_first (',' block_first)* ','? %optional_field;

    bed_entry = chrom_id first last (
                    name ( '\t' score ( '\t' strand ( '\t' thick_first (
                    '\t' thick_last ( '\t' item_rgb ( '\t' block_count (
                    '\t' block_sizes ( '\t' block_firsts )? )? )? )? )? )? )? )?
                    )? '\0';
    main := (bed_entry %finish_match)*;
}%%


%%write data;

# because brackets can fuck with ragel's parsing in unpredictable ways
typealias StringFieldVector Vector{StringField}
typealias NullableStringFieldVector Nullable{StringFieldVector}
typealias NullableStringField Nullable{StringField}

type BigBedDataParser
    state::Ragel.State

    # intermediate values used during parsing
    chrom_id::UInt32
    red::Float32
    green::Float32
    blue::Float32
    block_size_idx::Int
    block_first_idx::Int
    seq_names::Nullable{StringFieldVector}
    assumed_seqname::Nullable{StringField}

    function BigBedDataParser(input::BufferedInputStream;
                              seq_names::NullableStringFieldVector=NullableStringFieldVector(),
                              assumed_seqname::NullableStringField=NullableStringField())

        %% write init;

        return new(Ragel.State(cs, input), 0, 0.0, 0.0, 0.0, 1, 1, seq_names, assumed_seqname)
    end
end


Ragel.@generate_read_fuction("bigbed", BigBedDataParser, BEDInterval,
    begin
        %%write exec;
    end)


end # module BigBedDataParserImpl

