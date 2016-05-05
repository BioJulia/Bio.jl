# FASTA sequence types

immutable FASTA <: FileFormat end


"Metadata for FASTA sequence records containing just a `description` field"
type FASTAMetadata
    description::StringField
end


function FASTAMetadata()
    return FASTAMetadata(StringField())
end


function Base.(:(==))(a::FASTAMetadata, b::FASTAMetadata)
    return a.description == b.description
end


function Base.copy(metadata::FASTAMetadata)
    return FASTAMetadata(copy(metadata.description))
end


"FASTASeqRecord{S} is a `SeqRecord` for FASTA sequences of type `S`"
typealias FASTASeqRecord           SeqRecord{Sequence, FASTAMetadata}

"A `SeqRecord` type for FASTA DNA sequences"
typealias FASTADNASeqRecord       DNASeqRecord{FASTAMetadata}

"A `SeqRecord` type for FASTA RNA sequences"
typealias FASTARNASeqRecord       RNASeqRecord{FASTAMetadata}

"A `SeqRecord` type for FASTA amino acid sequences"
typealias FASTAAminoAcidSeqRecord AminoAcidSeqRecord{FASTAMetadata}

function Base.show{S}(io::IO, seqrec::SeqRecord{S, FASTAMetadata})
    write(io, ">", seqrec.name, " ", seqrec.metadata.description, "\n")
    show(io, seqrec.seq)
end

function Base.print{S}(io::IO, seqrec::SeqRecord{S,FASTAMetadata})
    write(io, ">", seqrec.name, " ", seqrec.metadata.description, "\n")
    print(io, seqrec.seq)
end


"Writes a FASTASeqRecord to an IO-stream (and obeys FASTAs max character constraint)"
function Base.write{T}(io::IO, seqrec::SeqRecord{T, FASTAMetadata})
    write(io, ">", seqrec.name)
    if !isempty(seqrec.metadata.description)
        write(io, " ", seqrec.metadata.description)
    end
    write(io, "\n")
    maxchars = 79
    counter = 1
    len = length(seqrec.seq)
    for nt in seqrec.seq
        show(io, nt)
        if counter % maxchars == 0 && counter < len
            write(io, "\n")
        end
        counter += 1
    end
    write(io, "\n")
end


%%{
    machine _fastaparser;

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
        yield = true;
        fbreak;
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
        %% write init;

        return new(Ragel.State(cs, input), BufferedOutputStream())
    end
end


function Base.eltype(::Type{FASTAParser})
    return FASTASeqRecord
end


function Base.open(input::BufferedInputStream, ::Type{FASTA})
    return FASTAParser(input)
end


typealias FASTAAnySeqRecord{S} SeqRecord{S, FASTAMetadata}

Ragel.@generate_read_fuction("_fastaparser", FASTAParser, FASTAAnySeqRecord,
    begin
        %% write exec;
    end)



