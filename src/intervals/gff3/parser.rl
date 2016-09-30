#

%%{
    machine gff3parser;

    action finish_match {
        Ragel.@anchor!
        Ragel.@yield ftargs
    }

    action count_line { state.linenum += 1 }
    action anchor { Ragel.@anchor! }

    action directive {
        Ragel.@copy_from_anchor!(input.directive)
        if startswith(input.directive, "gff-version")
            # ##gff-version 3.2.1
            input.version = VersionNumber(split(String(input.directive), r"\s+")[2])
        elseif startswith(input.directive, "sequence-region")
            # ##sequence-region seqid start end
            vals = split(String(input.directive), r"\s+")
            push!(input.sequence_regions,
                  Interval(vals[2], parse(Int, vals[3]), parse(Int, vals[4])))
        elseif startswith(input.directive, "FASTA")
            input.fasta_seen = true
            input.state.finished = true
            if input.entry_seen
                Ragel.@yield ftargs
            else
                fbreak;
            end
        end

        if input.save_directives
            if input.directive_count == length(input.directives)
                push!(input.directives, StringField())
            end
            copy!(input.directives[input.directive_count += 1], input.directive)
        end
    }

    action implicit_fasta {
        input.fasta_seen = true
        input.state.finished = true
        if input.entry_seen
            p -= 2
            Ragel.@yield ftargs
        else
            p -= 1
            fbreak;
        end
    }

    # interval
    action seqname {
        input.preceding_directives, input.directives =
            input.directives, input.preceding_directives
        input.preceding_directive_count = input.directive_count
        input.directive_count = 0;

        input.entry_seen = true
        empty!(output.metadata.attributes)
        Ragel.@copy_from_anchor!(output.seqname)
        if input.unescape_needed
            unescape!(output.metadata.source)
            input.unescape_needed = false
        end
    }
    action start   { output.first = Ragel.@int64_from_anchor! }
    action end     { output.last = Ragel.@int64_from_anchor! }
    action strand  { output.strand = convert(Strand, (Ragel.@char)) }

    # metadata
    action source     {
        Ragel.@copy_from_anchor!(output.metadata.source)
        if input.unescape_needed
            unescape!(output.metadata.source)
            input.unescape_needed = false
        end
    }
    action kind       { Ragel.@copy_from_anchor!(output.metadata.kind) }
    action score      { output.metadata.score = Nullable(Ragel.@float64_from_anchor!) }
    action nullscore  { output.metadata.score = Nullable{Float64}() }
    action phase      { output.metadata.phase = Nullable(Ragel.@int64_from_anchor!) }
    action nullphase  { output.metadata.phase = Nullable{Int}() }
    action attribute_key {
        Ragel.@copy_from_anchor!(input.key)
        if input.unescape_needed
            unescape!(output.metadata.source)
            input.unescape_needed = false
        end
    }
    action attribute_value {
        pushindex!(output.metadata.attributes, input.key,
                   input.state.stream.buffer, upanchor!(input.state.stream), p,
                   input.unescape_needed)
        input.unescape_needed = false
    }
    action check_percent {
        input.unescape_needed = input.unescape_needed || data[p] == UInt32('%')
    }

    newline        = '\r'? '\n' >count_line;
    hspace         = [ \t\v];
    blankline      = hspace* newline;
    # Just check that there's a digit. Actual validation happens when we parse
    # the float.
    floating_point = [ -~]* digit [ -~]*;

    comment   = '#' (any - newline - '#')* newline;
    directive = "##" ((any - newline)* >anchor %directive) newline;
    implicit_fasta = '>' >implicit_fasta;

    seqname    = ([a-zA-Z0-9.:^*$@!+_?\-|%] %check_percent)* >anchor %seqname;
    source     = ([ -~] %check_percent)* >anchor %source;
    kind       = [ -~]* >anchor %kind;
    start      = digit+ >anchor %start;
    end        = digit+ >anchor %end;
    score      = ((floating_point %score) | ('.' %nullscore)) >anchor;
    strand     = [+\-\.?] >strand;
    phase      = (([0-2] %phase) | ('.' %nullphase)) >anchor;

    attribute_char = ([ -~] - [=;,] %check_percent);
    attribute_key = attribute_char* >anchor %attribute_key;
    attribute_value = attribute_char* >anchor %attribute_value;
    attribute = attribute_key '=' attribute_value (',' attribute_value)*;
    attributes = (attribute ';')* attribute?;

    non_entry = blankline | directive | comment | implicit_fasta;

    gff3_entry = seqname '\t' source '\t' kind '\t' start '\t' end '\t'
                 score   '\t' strand '\t' phase '\t' attributes
                 newline non_entry*;

    main := non_entry* (gff3_entry %finish_match)*;
}%%

%% write data;

Ragel.@generate_read!_function(
    "gff3parser",
    GFF3Reader,
    GFF3Interval,
    begin
        %% write exec;
    end)
