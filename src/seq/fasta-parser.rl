%%{
    machine fastaparser;

    action finish_match {
        if seqtype(typeof(output)) == ReferenceSequence
            output.seq = ReferenceSequence(input.seqbuf.buffer, 1, length(input.seqbuf))
        elseif seqtype(typeof(output)) == BioSequence
            ET = predict(input.seqbuf.buffer, 1, length(input.seqbuf))
            if ET == typeof(output.seq)
                resize!(output.seq, length(input.seqbuf))
                encode_copy!(output.seq, 1, input.seqbuf.buffer, 1, length(input.seqbuf))
            else
                output.seq = ET(input.seqbuf.buffer, 1, length(input.seqbuf))
            end
        else
            resize!(output.seq, length(input.seqbuf))
            encode_copy!(output.seq, 1, input.seqbuf.buffer, 1, length(input.seqbuf))
        end
        empty!(input.seqbuf)
        Ragel.@yield ftargs
    }

    action count_line  { state.linenum += 1 }
    action mark        { Ragel.@anchor! }
    action identifier  { Ragel.@copy_from_anchor!(output.name) }
    action description { Ragel.@copy_from_anchor!(output.metadata.description) }
    action letters     { Ragel.@append_from_anchor!(input.seqbuf) }

    newline     = '\r'? '\n'     >count_line;
    hspace      = [ \t\v];
    whitespace  = space | newline;

    identifier  = (any - space)+            >mark  %identifier;
    description = ((any - hspace) [^\r\n]*) >mark  %description;
    letters     = (any - space - '>')+      >mark  %letters;
    sequence    = whitespace* letters? (whitespace+ letters)*;
    fasta_entry = '>' identifier (hspace+ description)? newline sequence whitespace*;

    main := whitespace* (fasta_entry %finish_match)**;
}%%

%% write data;

# FIXME: output type may be too loose
Ragel.@generate_read!_function(
    "fastaparser",
    FASTAParser,
    SeqRecord,
    begin
        %% write exec;
    end)
