%%{
    machine fastaparser;

    action finish_match {
        if seqtype(typeof(output)) == Sequence
            alphabet = predict(input.seqbuf.buffer, 1, length(input.seqbuf))
            ET = alphabet_type[alphabet]
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


"A type encapsulating the current state of a FASTA parser"
type FASTAParser <: AbstractParser
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}

    function FASTAParser(input::BufferedInputStream)
        return new(Ragel.State(fastaparser_start, input), BufferedOutputStream())
    end
end

function Base.eltype(::Type{FASTAParser})
    return FASTASeqRecord
end

function Base.eof(parser::FASTAParser)
    return eof(parser.state.stream)
end

function Base.open(input::BufferedInputStream, ::Type{FASTA})
    return FASTAParser(input)
end

Ragel.@generate_read_fuction(
    "fastaparser",
    FASTAParser,
    FASTASeqRecord,
    begin
        %% write exec;
    end)
