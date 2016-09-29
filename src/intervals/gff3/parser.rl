#

%%{
    machine gff3parser;

    action finish_match {
        Ragel.@anchor!
        Ragel.@yield ftargs
    }

    action count_line { input.state.linenum += 1 }
    action anchor { Ragel.@anchor! }

    action directive {
        directive = Ragel.@ascii_from_anchor!
        if startswith(directive, "##gff-version")
            # ##gff-version 3.2.1
            input.version = VersionNumber(split(directive, r"\s+")[2])
        elseif startswith(directive, "##sequence-region")
            # ##sequence-region seqid start end
            vals = split(directive, r"\s+")
            push!(input.sequence_regions,
                Interval(vals[2], parse(Int, vals[3]), parse(Int, vals[4])))
        else
            # TODO: record other directives
        end
    }

    # interval
    action seqname { Ragel.@copy_from_anchor!(output.seqname) }
    action start   { output.first = Ragel.@int64_from_anchor! }
    action end     { output.last = Ragel.@int64_from_anchor! }
    action strand  { output.strand = convert(Strand, (Ragel.@char)) }

    # metadata
    action source     { Ragel.@copy_from_anchor!(output.metadata.source) }
    action kind       { Ragel.@copy_from_anchor!(output.metadata.kind) }
    action score      { output.metadata.score = Nullable(Ragel.@float64_from_anchor!) }
    action nullscore  { output.metadata.score = Nullable{Float64}() }
    action phase      { output.metadata.phase = Nullable(Ragel.@int64_from_anchor!) }
    action nullphase  { output.metadata.phase = Nullable{Int}() }
    action attribute_key {
        Ragel.@copy_from_anchor!(input.key)
    }
    action attribute_value {
        pushindex!(output.metadata.attributes, input.key,
                   input.state.stream.buffer, upanchor!(input.state.stream), p)
    }

    newline        = '\r'? '\n' >count_line;
    hspace         = [ \t\v];
    blankline      = hspace* newline;
    # Just check that there's a digit. Actual validation happens when we parse
    # the float.
    floating_point = [ -~]* digit [ -~]*;

    comment   = '#' (any - newline - '#')* newline;
    directive = ("##" (any - newline)* newline) >anchor %directive;

    seqname    = [a-zA-Z0-9.:^*$@!+_?\-|%]* >anchor %seqname;
    source     = [ -~]* >anchor %source;
    kind       = [ -~]* >anchor %kind;
    start      = digit+ >anchor %start;
    end        = digit+ >anchor %end;
    score      = ((floating_point %score) | ('.' %nullscore)) >anchor;
    strand     = [+\-\.?] >strand;
    phase      = (([0-2] %phase) | ('.' %nullphase)) >anchor;

    attribute_char = [ -~] - [=;,];
    attribute_key = attribute_char* >anchor %attribute_key;
    attribute_value = attribute_char* >anchor %attribute_value;
    attribute = attribute_key '=' attribute_value (',' attribute_value)*;
    attributes = (attribute ';')* attribute?;

    gff3_entry = seqname '\t' source '\t' kind '\t' start '\t' end '\t'
                 score   '\t' strand '\t' phase '\t' attributes newline;

    main := (blankline | directive | comment | (gff3_entry %finish_match))*;
}%%

%% write data;

Ragel.@generate_read!_function(
    "gff3parser",
    GFF3Reader,
    GFF3Interval,
    begin
        %% write exec;
    end)
