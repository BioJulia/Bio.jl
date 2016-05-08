%%{
    machine fastqparser;

    action count_line {
        state.linenum += 1
    }

    action anchor { Ragel.@anchor! }

    action identifier   { Ragel.@copy_from_anchor!(output.name) }
    action description  { Ragel.@copy_from_anchor!(output.metadata.description) }
    action identifier2  { Ragel.@copy_from_anchor!(input.name2buf) }
    action description2 { Ragel.@copy_from_anchor!(input.desc2buf) }
    action letters      { Ragel.@append_from_anchor!(input.seqbuf) }
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
        encode_copy!(output.seq, 1, input.seqbuf.buffer, 1, input.seqbuf.position - 1)

        # quality
        encoding, input.quality_encodings =
            infer_quality_encoding(input.qualbuf.buffer, 1,
                                   input.qualbuf.position - 1,
                                   input.quality_encodings)

        decode_quality_string!(encoding, input.qualbuf.buffer,
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


"A type encapsulating the current state of a FASTQ parser"
type FASTQParser <: AbstractParser
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    qualbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    name2buf::StringField
    desc2buf::StringField
    qualcount::Int
    quality_encodings::QualityEncoding

    function FASTQParser(input::BufferedInputStream,
                         quality_encodings::QualityEncoding)
        return new(Ragel.State(fastqparser_start, input),
                   BufferedOutputStream(), BufferedOutputStream(),
                   StringField(), StringField(), 0, quality_encodings)
    end
end


function Base.eltype(::Type{FASTQParser})
    return FASTQSeqRecord
end

function Base.eof(parser::FASTQParser)
    return eof(parser.state.stream)
end


function Base.open(input::BufferedInputStream, ::Type{FASTQ};
                   quality_encodings::QualityEncoding=EMPTY_QUAL_ENCODING)
    return FASTQParser(input, quality_encodings)
end


Ragel.@generate_read_fuction("fastqparser", FASTQParser, FASTQSeqRecord,
    begin
        %% write exec;
    end)
