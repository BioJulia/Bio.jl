

# FASTQ sequence types

@doc """
Metadata for FASTQ sequence records containing a `description` field,
and a `quality` string corresponding to the sequence.

Quality scores are stored as integer Phred scores.
""" ->
type FASTQMetadata
    description::String
    quality::Vector{Int8}

    function FASTQMetadata(description, quality)
        return new(description, quality)
    end

    function FASTQ()
        return new("", Int8[])
    end
end


@doc """
A `SeqRecord` for FASTQ sequences.
""" ->
typealias FASTQSeqRecord DNASeqRecord{FASTQMetadata}


function Base.show(io::IO, seqrec::FASTQSeqRecord)
    write(io, "@", seqrec.name, " ", seqrec.metadata.description, "\n")
    for c in seqrec.seq
        show(io, c)
    end
    write(io, '\n')
    # print quality scores as a unicode bar chart
    for q in seqrec.metadata.quality
        if q <= 0
            write(io, '▁')
        elseif q <= 6
            write(io, '▂')
        elseif q <= 12
            write(io, '▃')
        elseif q <= 18
            write(io, '▄')
        elseif q <= 24
            write(io, '▅')
        elseif q <= 30
            write(io, '▆')
        elseif q <= 36
            write(io, '▇')
        else
            write(io, '█')
        end
    end
    write(io, '\n')
end


module FASTQParserImpl

import Bio.Seq: FASTQSeqRecord, QualityEncoding, EMPTY_QUAL_ENCODING
import Bio.Ragel
using Docile, Switch
export FASTQParser


%%{
    machine fastq;

    action yield {
        yield = true;
        fbreak;
    }

    action count_line {
        input.state.linenum += 1
    }

    action identifier_start {
        Ragel.@pushmark!
    }

    action identifier_end {
        firstpos = Ragel.@popmark!
        append!(input.namebuf, state.buffer, firstpos, p)
    }

    action description_start {
        Ragel.@pushmark!
    }

    action description_end {
        firstpos = Ragel.@popmark!
        append!(input.descbuf, state.buffer, firstpos, p)
    }

    action identifier2_start {
        Ragel.@pushmark!
    }

    action identifier2_end {
        firstpos = Ragel.@popmark!
        append!(input.name2buf, state.buffer, firstpos, p)
    }

    action description2_start {
        Ragel.@pushmark!
    }

    action description2_end {
        firstpos = Ragel.@popmark!
        append!(input.desc2buf, state.buffer, firstpos, p)
    }

    action letters_start {
        Ragel.@pushmark!
    }

    action letters_end {
        firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    }

    action qletters_start {
        Ragel.@pushmark!
    }

    action qletters_end {
        firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
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

    newline     = '\r'? '\n'     >count_line;
    hspace      = [ \t\v];
    whitespace  = newline | hspace;

    identifier  = (any - space)+     >identifier_start   %identifier_end;
    description = [^\r\n]+           >description_start  %description_end;

    identifier2  = (any - space)+    >identifier2_start  %identifier2_end;
    description2 = [^\r\n]+          >description2_start %description2_end;

    letters     = alpha+             >letters_start      %letters_end;
    sequence    = (newline+ letters)*;

    qletters    = ([!-~] when qlen_lt $inc_qual_count)+   >qletters_start %qletters_end;
    quality     = (newline+ qletters)*;

    fastq_entry = '@' when qlen_eq identifier ([ \t\v]+ description)?
                  sequence
                  newline+ '+' hspace* (identifier2 ([ \t\v]+ description2)?)?
                  quality newline+ whitespace*;

    main := whitespace* (fastq_entry %yield)*;
}%%


%% write data;


@doc """
A type encapsulating the current state of a FASTQ parser.
""" ->
type FASTQParser
    state::Ragel.State
    seqbuf::Ragel.Buffer{Uint8}
    qualbuf::Ragel.Buffer{Uint8}
    namebuf::Ragel.Buffer{Uint8}
    descbuf::Ragel.Buffer{Uint8}
    name2buf::Ragel.Buffer{Uint8}
    desc2buf::Ragel.Buffer{Uint8}
    qualcount::Int
    default_qual_encoding::QualityEncoding

    function FASTQParser(input::Union(IO, String, Vector{Uint8}),
                         default_qual_encoding=EMPTY_QUAL_ENCODING;
                         memory_map::Bool=false)
        %% write init;

        if memory_map
            if !isa(input, String)
                error("Parser must be given a file name in order to memory map.")
            end
            return new(Ragel.State(cs, input, true),
                       Ragel.Buffer{Uint8}(),
                       Ragel.Buffer{Uint8}(), Ragel.Buffer{Uint8}(),
                       Ragel.Buffer{Uint8}(), Ragel.Buffer{Uint8}(),
                       Ragel.Buffer{Uint8}(), 0, default_qual_encoding)
        else
            return new(Ragel.State(cs, input), Ragel.Buffer{Uint8}(),
                       Ragel.Buffer{Uint8}(), Ragel.Buffer{Uint8}(),
                       Ragel.Buffer{Uint8}(), Ragel.Buffer{Uint8}(),
                       Ragel.Buffer{Uint8}(), 0, default_qual_encoding)
        end
    end
end


function Ragel.ragelstate(parser::FASTQParser)
    return parser.state
end


function accept_state!(input::FASTQParser, output::FASTQSeqRecord)
    if length(input.seqbuf) != length(input.qualbuf)
        error("Error parsing FASTQ: sequence and quality scores must be of equal length")
    end
    output.name = input.namebuf
    output.metadata.description = input.descbuf
    output.seq = DNASequence(input.seqbuf.data, 1, input.seqbuf.pos - 1)

    encoding = infer_quality_encoding(input.qualbuf.data, 1,
                                      input.qualbuf.pos - 1,
                                      input.default_qual_encoding)
    input.default_qual_encoding = encoding
    output.metadata.quality = decode_quality_string(encoding, input.qualbuf.
                                                    1, input.qualbuf.pos - 1)

    input.namebuf = ""
    input.descbuf = ""
    empty!(input.seqbuf)
    empty!(input.qualbuf)
end


Ragel.@generate_read_fuction("fastq", FASTQParser, FASTQSeqRecord,
    begin
        @inbounds begin
            %% write exec;
        end
    end,
    begin
        accept_state!(input, output)
    end)

end # module FASTQParserImpl


using Bio.Seq.FASTQParserImpl

@doc """
An iterator over entries in a FASTQ file or stream.
""" ->
type FASTQIterator
    parser::FASTQParser

    # A type or function used to construct output sequence types
    default_qual_encoding::QualityEncoding
    isdone::Bool
    nextitem
end

@doc """
Parse a FASTQ file.

# Arguments
  * `filename::String`: Path of the FASTA file.
  * `qual_encoding::QualityEncoding`: assumed quality score encoding
    (Default: EMPTY_QUAL_ENCODING, i.e. no assumption)
  * `memory_map::Bool`: If true, attempt to memory map the file on supported
    platforms. (Default: `false`)

# Returns
An iterator over `SeqRecord`s contained in the file.
""" ->
function Base.read(filename::String, ::Type{FASTQ},
                   qual_encoding::QualityEncoding=EMPTY_QUAL_ENCODING;
                   memory_map=false)
    return FASTQIterator(FASTQParser(filename, memory_map=memory_map),
                         qual_encoding, false, nothing)
end

@doc """
Parse a FASTQ file.

# Arguments
  * `input::IO`: Input stream containing FASTQ data.
  * `qual_encoding::QualityEncoding`: assumed quality score encoding
    (Default: EMPTY_QUAL_ENCODING, i.e. no assumption)

# Returns
An iterator over `SeqRecord`s contained in the file.
""" ->
function Base.read(input::IO, ::Type{FASTQ},
                   qual_encoding::QualityEncoding=EMPTY_QUAL_ENCODING)
    return FASTQIterator(FASTQParser(input), qual_encoding, false, nothing)
end


function advance!(it::FASTQIterator)
    it.isdone = !FASTQParserImpl.advance!(it.parser)
    if !it.isdone
        if length(it.parser.seqbuf) != length(it.parser.qualbuf)
            error("Error parsing FASTQ: sequence and quality scores must be of equal length")
        end
        encoding = infer_quality_encoding(it.parser.qualbuf.data, 1,
                                          it.parser.qualbuf.pos - 1,
                                          it.default_qual_encoding)
        it.default_qual_encoding = encoding
        qscores = decode_quality_string(encoding, it.parser.qualbuf.data,
                                        1, it.parser.qualbuf.pos - 1)

        if (!isempty(it.parser.name2buf) && it.parser.namebuf != it.parser.name2buf) ||
           (!isempty(it.parser.desc2buf) && it.parser.descbuf != it.parser.desc2buf)
            error("Error parsing FASTQ: sequance and quality scores have non-matching identifiers")
        end

        it.nextitem =
            FASTQSeqRecord(takebuf_string(it.parser.namebuf),
                           DNASequence(it.parser.seqbuf.data, 1, it.parser.seqbuf.pos - 1),
                           FASTQMetadata(takebuf_string(it.parser.descbuf), qscores))
        empty!(it.parser.seqbuf)
        empty!(it.parser.qualbuf)
        empty!(it.parser.name2buf)
        empty!(it.parser.desc2buf)
    end
end


function start(it::FASTQIterator)
    advance!(it)
    return nothing
end


function next(it::FASTQIterator, state)
    item = it.nextitem
    advance!(it)
    return item, nothing
end


function done(it::FASTQIterator, state)
    return it.isdone
end


