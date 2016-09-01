%%{
    machine bedparser;

    action finish_match {
        input.block_size_idx = 1
        input.block_first_idx = 1
        Ragel.@anchor!
        Ragel.@yield ftargs
    }

    action count_line { input.state.linenum += 1 }
    action anchor { Ragel.@anchor! }
    action optional_field { output.metadata.used_fields += 1 }

    action seqname     { output.metadata.used_fields = 0; Ragel.@copy_from_anchor!(output.seqname) }
    action first       { output.first = 1 + Ragel.@int64_from_anchor! }
    action last        { output.last = Ragel.@int64_from_anchor! }
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

    newline      = '\r'? '\n'     >count_line;
    hspace       = [ \t\v];
    blankline    = hspace* newline;

    seqname      = [ -~]*   >anchor     %seqname;
    first        = digit+   >anchor     %first;
    last         = digit+   >anchor     %last;
    name         = [ -~]*   >anchor     %name  %optional_field;
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

    bed_entry = seqname '\t' first '\t' last (
                    '\t' name ( '\t' score ( '\t' strand ( '\t' thick_first (
                    '\t' thick_last ( '\t' item_rgb ( '\t' block_count (
                    '\t' block_sizes ( '\t' block_firsts )? )? )? )? )? )? )? )?
                    )?
                newline blankline*;

    main := blankline* (bed_entry %finish_match)*;
}%%

%% write data;

Ragel.@generate_read!_function(
    "bedparser",
    BEDReader,
    BEDInterval,
    begin
        %% write exec;
    end)
