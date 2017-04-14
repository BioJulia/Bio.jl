%%{
    machine fastqparser;

    action count_line {
        state.linenum += 1
    }

    action anchor { Ragel.@anchor!; }

    action identifier   { Ragel.@copy_from_anchor!(output.name); }
    action description  { Ragel.@copy_from_anchor!(output.metadata.description); }
    action identifier2  { Ragel.@copy_from_anchor!(input.name2buf); }
    action description2 { Ragel.@copy_from_anchor!(input.desc2buf); }
    action letters      { Ragel.@append_from_anchor!(input.seqbuf); }
    action qletters {
        Ragel.@append_from_anchor!(input.qualbuf)
        input.qualcount = 0
    }

    action qlen_lt {
        length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    }

    action qlen_eq {
        length(input.qualbuf) == length(input.seqbuf)
    }

    action inc_qual_count {
        input.qualcount += 1
    }

    action finish_match {
        if length(input.seqbuf) != length(input.qualbuf)
            error("Error parsing FASTQ: sequence and quality scores must be of equal length")
        end

        # if a name and s description are repeated (i.e. old-fashioned fastq)
        # make sure they match.
        if (!isempty(input.name2buf) && input.name2buf != output.name) ||
           (!isempty(input.desc2buf) && input.desc2buf != output.metadata.description)
            error("Error parsing FASTQ: sequence and quality scores have non-matching identifiers")
        end

        # sequence
        resize!(output.seq, input.seqbuf.position - 1)
        if !isnull(input.fill_ambiguous)
            byte = convert(UInt8, convert(Char, get(input.fill_ambiguous)))
            for i in 1:input.seqbuf.position-1
                b = input.seqbuf.buffer[i]
                if isambiguous(convert(DNA, convert(Char, b)))
                    input.seqbuf.buffer[i] = byte
                end
            end
        end
        encode_copy!(output.seq, 1, input.seqbuf.buffer, 1, input.seqbuf.position - 1)

        # quality
        check_quality_string(input.quality_encoding, input.qualbuf.buffer, 1, input.qualbuf.position - 1)
        decode_quality_string!(input.quality_encoding, input.qualbuf.buffer,
                               output.metadata.quality, 1,
                               input.qualbuf.position - 1)

        # reset temporaries for the next run
        empty!(input.qualbuf)
        empty!(input.seqbuf)
        empty!(input.name2buf)
        empty!(input.desc2buf)

        Ragel.@yield ftargs
    }

    newline     = '\r'? '\n'     >count_line;
    hspace      = [ \t\v];
    whitespace  = space | newline;

    identifier  = (any - space)+            >anchor  %identifier;
    description = ((any - space) [^\r\n]*)  >anchor  %description;

    identifier2  = (any - space)+           >anchor  %identifier2;
    description2 = ((any - space) [^\r\n]*) >anchor  %description2;

    letters     = [A-z]+                    >anchor  %letters;
    sequence    = letters? (newline+ letters)*;

    qletters    = ([!-~] when qlen_lt $inc_qual_count)+   >anchor %qletters;
    quality     = qletters? (newline+ qletters)*;

    fastq_entry = ('@' when qlen_eq) identifier (hspace+ description)?
                  newline sequence
                  newline+ '+' (identifier2 (hspace+ description2)?)?
                  newline quality newline+;

    main := whitespace* (fastq_entry %finish_match)**;
}%%

%% write data;

# FIXME: output type may be too loose
Ragel.@generate_read!_function(
    "fastqparser",
    FASTQReader,
    SeqRecord,
    begin
        %% write exec;
    end)
