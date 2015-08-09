

# FASTA sequence types

immutable FASTA <: FileFormat end


"Metadata for FASTA sequence records containing just a `description` field"
type FASTAMetadata
    description::String

    function FASTAMetadata(description)
        return new(description)
    end

    function FASTAMetadata()
        return new("")
    end
end


"FASTASeqRecord{S} is a `SeqRecord` for FASTA sequences of type `S`"
typealias FASTASeqRecord{S}       SeqRecord{S, FASTAMetadata}

"A `SeqRecord` type for FASTA DNA sequences"
typealias FASTADNASeqRecord       DNASeqRecord{FASTAMetadata}

"A `SeqRecord` type for FASTA RNA sequences"
typealias FASTARNASeqRecord       RNASeqRecord{FASTAMetadata}

"A `SeqRecord` type for FASTA amino acid sequences"
typealias FASTAAminoAcidSeqRecord AminoAcidSeqRecord{FASTAMetadata}


function Base.show(io::IO, seqrec::FASTASeqRecord)
    write(io, ">", seqrec.name, " ", seqrec.metadata.description, "\n")
    show(io, seqrec.seq)
end


"Writes a FASTASeqRecord to an IO-stream (and obeys FASTAs max character constraint)"
function Base.write(io::IO, seqrec::FASTASeqRecord)
    header = strip(string(">", seqrec.name, " ", seqrec.metadata.description))
    write(io, header, "\n")
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


module FASTAParserImpl

import Bio.Seq: FASTASeqRecord
import Bio.Ragel
using Switch
export FASTAParser


%%{
    machine fasta;

    action yield {
        yield = true;
        fbreak;
    }

    action count_line      { input.state.linenum += 1 }
    action pushmark        { Ragel.@pushmark! }
    action identifier  { input.namebuf = Ragel.@asciistring_from_mark!  }
    action description { input.descbuf = Ragel.@asciistring_from_mark! }
    action letters     { append!(input.seqbuf, state.buffer, (Ragel.@popmark!), p) }

    newline     = '\r'? '\n'     >count_line;
    hspace      = [ \t\v];
    whitespace  = space | newline;

    identifier  = (any - space)+            >pushmark  %identifier;
    description = ((any - hspace) [^\r\n]*) >pushmark  %description;
    letters     = (any - space - '>')+      >pushmark  %letters;
    sequence    = whitespace* letters? (whitespace+ letters)*;
    fasta_entry = '>' identifier (hspace+ description)? newline sequence whitespace*;

    main := whitespace* (fasta_entry %yield)*;
}%%


%% write data;


"A type encapsulating the current state of a FASTA parser"
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


"An iterator over entries in a FASTA file or stream."
type FASTAIterator
    parser::FASTAParser

    # A type or function used to construct output sequence types
    default_alphabet::Alphabet
    isdone::Bool
    nextitem
end

"""
Parse a FASTA file.

# Arguments
  * `filename::String`: Path of the FASTA file.
  * `alphabet::Alphabet`: Assumed alphabet for the sequences contained in the
      file. (Default: `DNA_ALPHABET`)
  * `memory_map::Bool`: If true, attempt to memory map the file on supported
    platforms. (Default: `false`)

# Returns
An iterator over `SeqRecord`s contained in the file.
"""
function Base.read(filename::String, ::Type{FASTA},
                   alphabet::Alphabet=DNA_ALPHABET; memory_map::Bool=false)
    return FASTAIterator(FASTAParser(filename, memory_map=memory_map),
                         alphabet, false, nothing)
end


"""
Parse a FASTA file.

# Arguments
  * `input::IO`: Input stream containing FASTA data.
  * `alphabet::Alphabet`: Assumed alphabet for the sequences contained in the
      file. (Default: DNA_ALPHABET)

# Returns
An iterator over `SeqRecord`s contained in the file.
"""
function Base.read(input::IO, ::Type{FASTA}, alphabet::Alphabet=DNA_ALPHABET)
    return FASTAIterator(FASTAParser(input), alphabet, false, nothing)
end


function Base.read(input::Cmd, ::Type{FASTA}, alphabet::Alphabet=DNA_ALPHABET)
    return FASTAIterator(FASTAParser(open(input, "r")[1]), alphabet, false, nothing)
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
