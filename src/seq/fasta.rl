

# FASTA sequence types

immutable FASTA <: FileFormat end


"Metadata for FASTA sequence records containing just a `description` field"
type FASTAMetadata
    description::StringField
end


function FASTAMetadata()
    return FASTAMetadata(StringField())
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


module FASTAParserImpl

import Bio.Ragel
using Bio: AbstractParser, StringField
import Bio.Seq
using Bio.Seq: Alphabet, DNA_ALPHABET, infer_alphabet, alphabet_type,
               Sequence, FASTA, FASTAMetadata, SeqRecord, FASTASeqRecord, seqtype
using BufferedStreams
using Switch


%%{
    machine fasta;

    action finish_match {
        if seqtype(typeof(output)) == Sequence
            alphabet = infer_alphabet(input.seqbuf.buffer, 1,
                                       length(input.seqbuf), input.default_alphabet)
            ET = alphabet_type[alphabet]
            if ET == typeof(output.seq)
                copy!(output.seq, input.seqbuf.buffer, 1, length(input.seqbuf))
            else
                output.seq = ET(input.seqbuf.buffer, 1, length(input.seqbuf),
                                mutable=true)
            end
            input.default_alphabet = alphabet
        else
            copy!(output.seq, input.seqbuf.buffer, 1, length(input.seqbuf))
        end
        empty!(input.seqbuf)
        yield = true;
        fbreak;
    }

    action count_line  { state.linenum += 1 }
    action mark        { Ragel.anchor!(state, p) }
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

    main := whitespace* (fasta_entry %finish_match)*;
}%%


%% write data;


"A type encapsulating the current state of a FASTA parser"
type FASTAParser <: AbstractParser
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStreamSource}
    default_alphabet::Alphabet

    function FASTAParser(input::BufferedInputStream)
        %% write init;

        return new(Ragel.State(cs, input), BufferedOutputStream(), DNA_ALPHABET)
    end
end


function Base.eltype(::Type{FASTAParser})
    return FASTASeqRecord
end


function Base.open(input::BufferedInputStream, ::Type{FASTA})
    return FASTAParser(input)
end


typealias FASTAAnySeqRecord{S} SeqRecord{S, FASTAMetadata}

Ragel.@generate_read_fuction("fasta", FASTAParser, FASTAAnySeqRecord,
    begin
        %% write exec;
    end)


end # module FASTAParserImpl



