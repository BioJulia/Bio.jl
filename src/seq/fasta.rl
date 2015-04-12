

# FASTA sequence types

@doc """
Metadata for FASTA sequence records containing just a `description` field.
""" ->
type FASTAMetadata
    description::String

    function FASTAMetadata(description)
        return new(description)
    end

    function FASTAMetadata()
        return new("")
    end
end


@doc """
FASTASeqRecord{S} is a `SeqRecord` for FASTA sequences of type `S`.
""" ->
typealias FASTASeqRecord{S}       SeqRecord{S, FASTAMetadata}

@doc """
A `SeqRecord` type for FASTA DNA sequences.
""" ->
typealias FASTADNASeqRecord       DNASeqRecord{FASTAMetadata}

@doc """
A `SeqRecord` type for FASTA RNA sequences.
""" ->
typealias FASTARNASeqRecord       RNASeqRecord{FASTAMetadata}

@doc """
A `SeqRecord` type for FASTA amino acid sequences.
""" ->
typealias FASTAAminoAcidSeqRecord AminoAcidSeqRecord{FASTAMetadata}


function Base.show(io::IO, seqrec::FASTASeqRecord)
    write(io, ">", seqrec.name, " ", seqrec.metadata.description, "\n")
    show(io, seqrec.seq)
end


module FASTAParserImpl

import Bio.Seq: FASTASeqRecord
import Bio.Ragel
using Docile, Switch
export FASTAParser


%%{
    machine fasta;

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
        input.namebuf = bytestring(Ragel.@spanfrom firstpos)
    }

    action description_start {
        Ragel.@pushmark!
    }

    action description_end {
        firstpos = Ragel.@popmark!
        input.descbuf = bytestring(Ragel.@spanfrom firstpos)
    }

    action letters_start {
        Ragel.@pushmark!
    }

    action letters_end {
        firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    }


    newline     = '\r'? '\n'     >count_line;
    hspace      = [ \t\v];
    whitespace  = newline | hspace;

    identifier  = (any - space)+ >identifier_start  %identifier_end;
    description = [^\r\n]+       >description_start %description_end;
    letters     = alpha+         >letters_start     %letters_end;
    sequence    = whitespace* letters? (newline+ whitespace* letters (hspace+ letters)*)*;
    fasta_entry = '>' identifier ( hspace+ description )? newline sequence whitespace*;

    main := whitespace* (fasta_entry %yield)*;
}%%


%% write data;


@doc """
A type encapsulating the current state of a FASTA parser.
""" ->
type FASTAParser
    state::Ragel.State
    seqbuf::Ragel.Buffer
    namebuf::String
    descbuf::String

    function FASTAParser(input::Union(IO, String, Vector{Uint8});
                         memory_map::Bool=false)
        %% write init;

        if memory_map
            if !isa(input, String)
                error("Parser must be given a file name in order to memory map.")
            end
            return new(Ragel.State(cs, input, true),
                       Ragel.Buffer{Uint8}(), "", "")
        else
            return new(Ragel.State(cs, input), Ragel.Buffer{Uint8}(), "", "")
        end
    end
end


function Ragel.ragelstate(parser::FASTAParser)
    return parser.state
end


function accept_state!{S}(input::FASTAParser, output::FASTASeqRecord{S})
    output.name = input.namebuf
    output.metadata.description = input.descbuf
    output.seq = S(input.seqbuf.data, 1, input.seqbuf.pos - 1)

    input.namebuf = ""
    input.descbuf = ""
    empty!(input.seqbuf)
end


Ragel.@generate_read_fuction("fasta", FASTAParser, FASTASeqRecord,
    begin
        @inbounds begin
            %% write exec;
        end
    end,
    begin
        accept_state!(input, output)
    end)

end # module FASTAParserImpl


using Bio.Seq.FASTAParserImpl


@doc """
An iterator over entries in a FASTA file or stream.
""" ->
type FASTAIterator
    parser::FASTAParser

    # A type or function used to construct output sequence types
    default_alphabet::Alphabet
    isdone::Bool
    nextitem
end

@doc """
Parse a FASTA file.

# Arguments
  * `filename::String`: Path of the FASTA file.
  * `alphabet::Alphabet`: Assumed alphabet for the sequences contained in the
      file. (Default: `DNA_ALPHABET`)
  * `memory_map::Bool`: If true, attempt to memory map the file on supported
    platforms. (Default: `false`)

# Returns
An iterator over `SeqRecord`s contained in the file.
""" ->
function Base.read(filename::String, ::Type{FASTA},
                   alphabet::Alphabet=DNA_ALPHABET; memory_map::Bool=false)
    return FASTAIterator(FASTAParser(filename, memory_map=memory_map),
                         alphabet, false, nothing)
end


@doc """
Parse a FASTA file.

# Arguments
  * `input::IO`: Input stream containing FASTA data.
  * `alphabet::Alphabet`: Assumed alphabet for the sequences contained in the
      file. (Default: DNA_ALPHABET)

# Returns
An iterator over `SeqRecord`s contained in the file.
""" ->
function Base.read(input::IO, ::Type{FASTA}, alphabet::Alphabet=DNA_ALPHABET)
    return FASTAIterator(FASTAParser(input), alphabet, false, nothing)
end


function advance!(it::FASTAIterator)
    it.isdone = !FASTAParserImpl.advance!(it.parser)
    if !it.isdone
        alphabet = infer_alphabet(it.parser.seqbuf.data, 1, it.parser.seqbuf.pos - 1,
                                  it.default_alphabet)
        S = alphabet_type[alphabet]
        it.default_alphabet = alphabet
        it.nextitem =
            FASTASeqRecord{S}(it.parser.namebuf,
                              S(it.parser.seqbuf.data, 1, it.parser.seqbuf.pos - 1, true),
                              FASTAMetadata(it.parser.descbuf))
        empty!(it.parser.seqbuf)
    end
end


function start(it::FASTAIterator)
    advance!(it)
    return nothing
end


function next(it::FASTAIterator, state)
    item = it.nextitem
    advance!(it)
    return item, nothing
end


function done(it::FASTAIterator, state)
    return it.isdone
end

