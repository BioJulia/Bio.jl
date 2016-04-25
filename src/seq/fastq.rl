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


function Base.(:(==))(a::FASTQMetadata, b::FASTQMetadata)
    return a.description == b.description && a.quality == b.quality
end


function Base.copy(metadata::FASTQMetadata)
    return FASTQMetadata(copy(metadata.description), copy(metadata.quality))
end


"A `SeqRecord` for FASTQ sequences"
typealias FASTQSeqRecord DNASeqRecord{FASTQMetadata}


function FASTQSeqRecord()
    return FASTQSeqRecord(StringField(), DNASequence(), FASTQMetadata())
end

"""
Trim a FASTQSeqRecord in-place, including the quality score.
"""
function trim!(seqrec::FASTQSeqRecord, i::UnitRange)
    seqrec.metadata.quality = seqrec.metadata.quality[i]
    seqrec.seq = seqrec.seq[i]
end

"""
Slice a FASTQSeqRecord, including the quality score.
"""
function getindex(seqrec::FASTQSeqRecord, i::UnitRange)
    metadata = copy(seqrec.metadata)
    metadata.quality = metadata.quality[i]
    return FASTQSeqRecord(seqrec.name, seqrec.seq[i], metadata)
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

skipempty: If true, and the sequence is empty, then `write` becomes a no-op.
emptytoN: If true, and the sequence is empty, a single N is printed in place of
          the empty sequence to preserve FASTQ file validity.
"""
function Base.write(io::IO, seqrec::FASTQSeqRecord; offset::Integer=-1,
                    qualheader::Bool=false, skipempty::Bool=false,
                    emptytoN::Bool=false)

    # Skip empty records if we have been told to do so.
    if skipempty && length(seqrec) == 0
        return
    end

    # choose offset automatically
    if offset < 0
        if !isempty(seqrec.metadata.quality) && minimum(seqrec.metadata.quality) < 0
            offset = 64 # solexa quality offset
        else
            offset = 33  # sanger
        end
    end

    write(io, "@", seqrec.name)
    if !isempty(seqrec.metadata.description)
        write(io, " ", seqrec.metadata.description)
    end
    write(io, "\n")

    if emptytoN && length(seqrec) == 0
        write(io, "N")
    end
    for c in seqrec.seq
        show(io, c)
    end
    write(io, "\n")

    write(io, "+")
    if qualheader
        write(io, seqrec.name)
        if !isempty(seqrec.metadata.description)
            write(io, " ", seqrec.metadata.description)
        end
    end
    write(io, "\n")

    if emptytoN && length(seqrec) == 0
        # Write the char the corresponds to PHRED score of 40 in current
        # encoding.
        write(io, Char(offset + 40))
    end
    for q in seqrec.metadata.quality
        write(io, Char(q + offset))
    end
    write(io, "\n")
end


%%{
    machine _fastqparser;

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
        yield = true;
        fbreak;
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
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStreamSource}
    qualbuf::BufferedOutputStream{BufferedStreams.EmptyStreamSource}
    name2buf::StringField
    desc2buf::StringField
    qualcount::Int
    quality_encodings::QualityEncoding

    function FASTQParser(input::BufferedInputStream,
                         quality_encodings::QualityEncoding)
        %% write init;

        return new(Ragel.State(cs, input),
                   BufferedOutputStream(), BufferedOutputStream(),
                   StringField(), StringField(), 0, quality_encodings)
    end
end


function Base.eltype(::Type{FASTQParser})
    return FASTQSeqRecord
end


function Base.open(input::BufferedInputStream, ::Type{FASTQ};
                   quality_encodings::QualityEncoding=EMPTY_QUAL_ENCODING)
    return FASTQParser(input, quality_encodings)
end


Ragel.@generate_read_fuction("_fastqparser", FASTQParser, FASTQSeqRecord,
    begin
        %% write exec;
    end)



