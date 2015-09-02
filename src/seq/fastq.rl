

# FASTQ sequence types

immutable FASTQ <: FileFormat end


"""
Metadata for FASTQ sequence records containing a `description` field,
and a `quality` string corresponding to the sequence.

Quality scores are stored as integer Phred scores.
"""
type FASTQMetadata
    description::StringField
    quality::Vector{Int8}
end


function FASTQMetadata()
    return FASTQMetadata(StringField(), Int8[])
end


"A `SeqRecord` for FASTQ sequences"
typealias FASTQSeqRecord DNASeqRecord{FASTQMetadata}


function FASTQSeqRecord()
    return FASTQSeqRecord(StringField(), DNASequence(mutable=true),
                          FASTQMetadata())
end


"""
Show a `FASTQSeqRecord` to `io`, with graphical display of quality scores.
"""
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

"""
Write a `FASTQSeqRecord` to `io`, as a valid FASTQ record.
"""
function Base.write(io::IO, seqrec::FASTQSeqRecord; offset::Integer=33,
                    qualheader::Bool=false)
    write(io, "@", seqrec.name, " ", seqrec.metadata.description, "\n")

    for c in seqrec.seq
        show(io, c)
    end
    write(io, "\n")

    if qualheader
        write(io, "+", seqrec.name, " ", seqrec.metadata.description, "\n")
    else
        write(io, "+\n")
    end

    for q in seqrec.metadata.quality
        write(io, char(q + offset))
    end
    write(io, "\n")
end


module FASTQParserImpl

using Bio: StringField
using Bio.Seq: FASTQSeqRecord, QualityEncoding, EMPTY_QUAL_ENCODING,
               infer_quality_encoding, decode_quality_string!
import Bio.Ragel
using BufferedStreams
using Switch
export FASTQParser


%%{
    machine fastq;

    action count_line {
        state.linenum += 1
    }

    action anchor { Ragel.anchor!(state, p) }

    action identifier   { copy!(output.name, state.stream.buffer, Ragel.upanchor!(state), p) }
    action description  { copy!(output.metadata.description, state.stream.buffer, Ragel.upanchor!(state), p) }
    action identifier2  { append!(input.name2buf, state.stream.buffer, Ragel.upanchor!(state), p) }
    action description2 { append!(input.desc2buf, state.stream.buffer, Ragel.upanchor!(state), p) }
    action letters { append!(input.seqbuf, state.stream.buffer, Ragel.upanchor!(state), p) }
    action qletters {
        append!(input.qualbuf, state.stream.buffer, Ragel.upanchor!(state), p)
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

    action begin_match {
    }

    action finish_match {
        if length(input.seqbuf) != length(input.qualbuf)
            error("Error parsing FASTQ: sequence and quality scores must be of equal length")
        end

        # TODO: test name2buf and desc2buf

        # sequence
        copy!(output.seq, input.seqbuf.buffer, 1, input.seqbuf.position - 1)

        # quality
        encoding = infer_quality_encoding(input.qualbuf.buffer, 1,
                                          input.qualbuf.position - 1,
                                          input.default_qual_encoding)
        input.default_qual_encoding = encoding
        decode_quality_string!(encoding, input.qualbuf.buffer,
                               output.metadata.quality, 1,
                               input.qualbuf.position - 1)

        empty!(input.qualbuf)
        empty!(input.seqbuf)
        empty!(input.name2buf)
        empty!(input.desc2buf)
        yield = true;
        fbreak;
    }

    newline     = '\r'? '\n'     >count_line;
    hspace      = [ \t\v];
    whitespace  = space | newline;

    identifier  = (any - space)+           >anchor  %identifier;
    description = ((any - space) [^\r\n]*) >anchor  %description;

    identifier2  = (any - space)+           >anchor  %identifier2;
    description2 = ((any - space) [^\r\n]*) >anchor  %description2;

    letters     = [A-z]+                   >anchor  %letters;
    sequence    = letters? (newline+ letters)*;

    qletters    = ([!-~] when qlen_lt $inc_qual_count)+   >anchor %qletters;
    quality     = qletters? (newline+ qletters)*;

    fastq_entry = ('@' when qlen_eq) identifier (hspace+ description)?
                  newline sequence
                  newline+ '+' (identifier2 (hspace+ description2)?)?
                  newline quality newline+;

    main := whitespace* (fastq_entry >begin_match %finish_match)*;
}%%


%% write data;


"A type encapsulating the current state of a FASTQ parser"
type FASTQParser
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStreamSource}
    qualbuf::BufferedOutputStream{BufferedStreams.EmptyStreamSource}
    name2buf::StringField
    desc2buf::StringField
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
                       BufferedOutputStream(), BufferedOutputStream(),
                       StringField(), StringField(), 0,
                       default_qual_encoding)
        else
            return new(Ragel.State(cs, input), BufferedOutputStream(),
                       BufferedOutputStream(), StringField(), StringField(),
                       0, default_qual_encoding)
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
    output.seq = DNASequence(input.seqbuf.buffer, 1, input.seqbuf.position - 1)

    encoding = infer_quality_encoding(input.qualbuf.buffer, 1,
                                      input.qualbuf.position - 1,
                                      input.default_qual_encoding)
    input.default_qual_encoding = encoding
    output.metadata.quality = decode_quality_string(encoding, input.qualbuf.
                                                    1, input.qualbuf.position - 1)

    empty!(input.seqbuf)
    empty!(input.qualbuf)
end


Ragel.@generate_read_fuction("fastq", FASTQParser, FASTQSeqRecord,
    begin
        %% write exec;
    end)

end # module FASTQParserImpl


using Bio.Seq.FASTQParserImpl

"An iterator over entries in a FASTQ file or stream"
type FASTQIterator
    parser::FASTQParser

    # A type or function used to construct output sequence types
    default_qual_encoding::QualityEncoding
    isdone::Bool
    nextitem
end

"""
Parse a FASTQ file.

# Arguments
  * `filename::String`: Path of the FASTA file.
  * `qual_encoding::QualityEncoding`: assumed quality score encoding
    (Default: EMPTY_QUAL_ENCODING, i.e. no assumption)
  * `memory_map::Bool`: If true, attempt to memory map the file on supported
    platforms. (Default: `false`)

# Returns
An iterator over `SeqRecord`s contained in the file.
"""
function Base.read(filename::String, ::Type{FASTQ},
                   qual_encoding::QualityEncoding=EMPTY_QUAL_ENCODING;
                   memory_map=false)
    return FASTQIterator(FASTQParser(filename, memory_map=memory_map),
                         qual_encoding, false, nothing)
end

"""
Parse a FASTQ file.

# Arguments
  * `input::IO`: Input stream containing FASTQ data.
  * `qual_encoding::QualityEncoding`: assumed quality score encoding
    (Default: EMPTY_QUAL_ENCODING, i.e. no assumption)

# Returns
An iterator over `SeqRecord`s contained in the file.
"""
function Base.read(input::IO, ::Type{FASTQ},
                   qual_encoding::QualityEncoding=EMPTY_QUAL_ENCODING)
    return FASTQIterator(FASTQParser(input), qual_encoding, false, nothing)
end


function Base.read(input::Cmd, ::Type{FASTQ},
                   qual_encoding::QualityEncoding=EMPTY_QUAL_ENCODING)
    return FASTQIterator(FASTQParser(open(input, "r")[1]), qual_encoding, false, nothing)
end


function advance!(it::FASTQIterator)
    it.isdone = !FASTQParserImpl.advance!(it.parser)
    if !it.isdone
        if length(it.parser.seqbuf) != length(it.parser.qualbuf)
            error("Error parsing FASTQ: sequence and quality scores must be of equal length")
        end
        encoding = infer_quality_encoding(it.parser.qualbuf.buffer, 1,
                                          it.parser.qualbuf.position - 1,
                                          it.default_qual_encoding)
        it.default_qual_encoding = encoding
        qscores = decode_quality_string(encoding, it.parser.qualbuf.buffer,
                                        1, it.parser.qualbuf.position - 1)

        if (!isempty(it.parser.name2buf) && it.parser.name2buf != it.parser.namebuf) ||
           (!isempty(it.parser.desc2buf) && it.parser.desc2buf != it.parser.descbuf)
            error("Error parsing FASTQ: sequence and quality scores have non-matching identifiers")
        end

        seq = DNASequence(it.parser.seqbuf.buffer, 1, it.parser.seqbuf.position - 1)
        it.nextitem =
            FASTQSeqRecord(it.parser.namebuf, seq,
                           FASTQMetadata(it.parser.descbuf, qscores))

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
