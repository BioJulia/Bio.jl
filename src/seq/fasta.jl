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


const fasta_start  = convert(Int , 10)
const fasta_first_final  = convert(Int , 10)
const fasta_error  = convert(Int , 0)
const fasta_en_main  = convert(Int , 10)
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
        begin
cs = fasta_start;
	end
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
            begin
if p == pe 
	@goto _test_eof

end
@switch cs  begin
    @case 10
@goto st_case_10
@case 0
@goto st_case_0
@case 1
@goto st_case_1
@case 2
@goto st_case_2
@case 3
@goto st_case_3
@case 4
@goto st_case_4
@case 5
@goto st_case_5
@case 11
@goto st_case_11
@case 12
@goto st_case_12
@case 6
@goto st_case_6
@case 13
@goto st_case_13
@case 14
@goto st_case_14
@case 15
@goto st_case_15
@case 7
@goto st_case_7
@case 16
@goto st_case_16
@case 8
@goto st_case_8
@case 9
@goto st_case_9

end
@goto st_out
@label ctr0
begin
	input.state.linenum += 1
    
end
@goto st10
@label st10
p+= 1;
	if p == pe 
	@goto _test_eof10

end
@label st_case_10
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st10
end
@case 10
begin
@goto ctr0
end
@case 11
begin
@goto st10
end
@case 13
begin
@goto st1
end
@case 32
begin
@goto st10
end
@case 62
begin
@goto st2
end

end
begin
@goto st0
end
@label st_case_0
@label st0
cs = 0;
	@goto _out
@label st1
p+= 1;
	if p == pe 
	@goto _test_eof1

end
@label st_case_1
if ( data[1 + p ]) == 10 
	begin
@goto ctr0
end

end
begin
@goto st0
end
@label ctr20
begin
	yield = true;
        begin
	p+= 1; cs = 2; @goto _out

end

    
end
@goto st2
@label ctr27
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
begin
	yield = true;
        begin
	p+= 1; cs = 2; @goto _out

end

    
end
@goto st2
@label st2
p+= 1;
	if p == pe 
	@goto _test_eof2

end
@label st_case_2
if ( data[1 + p ]) == 32 
	begin
@goto st0
end

end
if ( data[1 + p ]) < 14 
	begin
if 9 <= ( data[1 + p ]) 
	begin
@goto st0
end

end
end

elseif ( ( data[1 + p ]) > 31  )
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr2
end

end
end

else
	begin
@goto ctr2
end

end
begin
@goto ctr2
end
@label ctr2
begin
	Ragel.@pushmark!
    
end
@goto st3
@label st3
p+= 1;
	if p == pe 
	@goto _test_eof3

end
@label st_case_3
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr4
end
@case 10
begin
@goto ctr5
end
@case 11
begin
@goto ctr4
end
@case 12
begin
@goto st0
end
@case 13
begin
@goto ctr6
end
@case 32
begin
@goto ctr4
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto st3
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto st3
end

end
begin
@goto st3
end
@label ctr4
begin
	firstpos = Ragel.@popmark!
        input.namebuf = bytestring(Ragel.@spanfrom firstpos)
    
end
@goto st4
@label st4
p+= 1;
	if p == pe 
	@goto _test_eof4

end
@label st_case_4
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr8
end
@case 11
begin
@goto ctr8
end
@case 12
begin
@goto ctr7
end
@case 32
begin
@goto ctr8
end

end
if ( data[1 + p ]) < 14 
	begin
if 10 <= ( data[1 + p ]) 
	begin
@goto st0
end

end
end

elseif ( ( data[1 + p ]) > 31  )
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr7
end

end
end

else
	begin
@goto ctr7
end

end
begin
@goto ctr7
end
@label ctr7
begin
	Ragel.@pushmark!
    
end
@goto st5
@label st5
p+= 1;
	if p == pe 
	@goto _test_eof5

end
@label st_case_5
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr10
end
@case 13
begin
@goto ctr11
end

end
if ( data[1 + p ]) > 12 
	begin
if 14 <= ( data[1 + p ]) 
	begin
@goto st5
end

end
end

elseif ( ( data[1 + p ]) >= 11  )
	begin
@goto st5
end

end
begin
@goto st5
end
@label ctr14
begin
	input.state.linenum += 1
    
end
@goto st11
@label ctr5
begin
	firstpos = Ragel.@popmark!
        input.namebuf = bytestring(Ragel.@spanfrom firstpos)
    
end
begin
	input.state.linenum += 1
    
end
@goto st11
@label ctr10
begin
	firstpos = Ragel.@popmark!
        input.descbuf = bytestring(Ragel.@spanfrom firstpos)
    
end
begin
	input.state.linenum += 1
    
end
@goto st11
@label st11
p+= 1;
	if p == pe 
	@goto _test_eof11

end
@label st_case_11
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st11
end
@case 10
begin
@goto ctr12
end
@case 11
begin
@goto st11
end
@case 13
begin
@goto st6
end
@case 32
begin
@goto st11
end
@case 62
begin
@goto ctr20
end

end
if ( data[1 + p ]) > 90 
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto ctr21
end

end
end

elseif ( ( data[1 + p ]) >= 65  )
	begin
@goto ctr21
end

end
begin
@goto st0
end
@label ctr12
begin
	input.state.linenum += 1
    
end
@goto st12
@label ctr25
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
@goto st12
@label st12
p+= 1;
	if p == pe 
	@goto _test_eof12

end
@label st_case_12
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st12
end
@case 10
begin
@goto ctr12
end
@case 11
begin
@goto st12
end
@case 13
begin
@goto st6
end
@case 32
begin
@goto st12
end
@case 62
begin
@goto ctr20
end

end
if ( data[1 + p ]) > 90 
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto ctr23
end

end
end

elseif ( ( data[1 + p ]) >= 65  )
	begin
@goto ctr23
end

end
begin
@goto st0
end
@label ctr26
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
@goto st6
@label st6
p+= 1;
	if p == pe 
	@goto _test_eof6

end
@label st_case_6
if ( data[1 + p ]) == 10 
	begin
@goto ctr12
end

end
begin
@goto st0
end
@label ctr23
begin
	Ragel.@pushmark!
    
end
@goto st13
@label st13
p+= 1;
	if p == pe 
	@goto _test_eof13

end
@label st_case_13
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr24
end
@case 10
begin
@goto ctr25
end
@case 11
begin
@goto ctr24
end
@case 13
begin
@goto ctr26
end
@case 32
begin
@goto ctr24
end
@case 62
begin
@goto ctr27
end

end
if ( data[1 + p ]) > 90 
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto st13
end

end
end

elseif ( ( data[1 + p ]) >= 65  )
	begin
@goto st13
end

end
begin
@goto st0
end
@label ctr24
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
@goto st14
@label st14
p+= 1;
	if p == pe 
	@goto _test_eof14

end
@label st_case_14
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st14
end
@case 10
begin
@goto ctr13
end
@case 11
begin
@goto st14
end
@case 13
begin
@goto st7
end
@case 32
begin
@goto st14
end
@case 62
begin
@goto ctr20
end

end
if ( data[1 + p ]) > 90 
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto ctr23
end

end
end

elseif ( ( data[1 + p ]) >= 65  )
	begin
@goto ctr23
end

end
begin
@goto st0
end
@label ctr13
begin
	input.state.linenum += 1
    
end
@goto st15
@label ctr32
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
@goto st15
@label st15
p+= 1;
	if p == pe 
	@goto _test_eof15

end
@label st_case_15
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st15
end
@case 10
begin
@goto ctr13
end
@case 11
begin
@goto st15
end
@case 13
begin
@goto st7
end
@case 32
begin
@goto st15
end
@case 62
begin
@goto ctr20
end

end
begin
@goto st0
end
@label st7
p+= 1;
	if p == pe 
	@goto _test_eof7

end
@label st_case_7
if ( data[1 + p ]) == 10 
	begin
@goto ctr13
end

end
begin
@goto st0
end
@label ctr21
begin
	Ragel.@pushmark!
    
end
@goto st16
@label st16
p+= 1;
	if p == pe 
	@goto _test_eof16

end
@label st_case_16
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr32
end
@case 10
begin
@goto ctr25
end
@case 11
begin
@goto ctr32
end
@case 13
begin
@goto ctr26
end
@case 32
begin
@goto ctr32
end
@case 62
begin
@goto ctr27
end

end
if ( data[1 + p ]) > 90 
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto st16
end

end
end

elseif ( ( data[1 + p ]) >= 65  )
	begin
@goto st16
end

end
begin
@goto st0
end
@label ctr6
begin
	firstpos = Ragel.@popmark!
        input.namebuf = bytestring(Ragel.@spanfrom firstpos)
    
end
@goto st8
@label ctr11
begin
	firstpos = Ragel.@popmark!
        input.descbuf = bytestring(Ragel.@spanfrom firstpos)
    
end
@goto st8
@label st8
p+= 1;
	if p == pe 
	@goto _test_eof8

end
@label st_case_8
if ( data[1 + p ]) == 10 
	begin
@goto ctr14
end

end
begin
@goto st0
end
@label ctr8
begin
	Ragel.@pushmark!
    
end
@goto st9
@label st9
p+= 1;
	if p == pe 
	@goto _test_eof9

end
@label st_case_9
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr8
end
@case 10
begin
@goto ctr10
end
@case 11
begin
@goto ctr8
end
@case 12
begin
@goto ctr7
end
@case 13
begin
@goto ctr11
end
@case 32
begin
@goto ctr8
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr7
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto ctr7
end

end
begin
@goto ctr7
end
@label st_out
@label _test_eof10
cs = 10; @goto _test_eof
@label _test_eof1
cs = 1; @goto _test_eof
@label _test_eof2
cs = 2; @goto _test_eof
@label _test_eof3
cs = 3; @goto _test_eof
@label _test_eof4
cs = 4; @goto _test_eof
@label _test_eof5
cs = 5; @goto _test_eof
@label _test_eof11
cs = 11; @goto _test_eof
@label _test_eof12
cs = 12; @goto _test_eof
@label _test_eof6
cs = 6; @goto _test_eof
@label _test_eof13
cs = 13; @goto _test_eof
@label _test_eof14
cs = 14; @goto _test_eof
@label _test_eof15
cs = 15; @goto _test_eof
@label _test_eof7
cs = 7; @goto _test_eof
@label _test_eof16
cs = 16; @goto _test_eof
@label _test_eof8
cs = 8; @goto _test_eof
@label _test_eof9
cs = 9; @goto _test_eof
@label _test_eof
begin
end
if p == eof 
	begin
@switch cs  begin
    @case 11
@case 12
@case 14
@case 15
begin
	yield = true;
        begin
	p+= 1; cs = 0; @goto _out

end

    
end

	break;
	@case 13
@case 16
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
begin
	yield = true;
        begin
	p+= 1; cs = 0; @goto _out

end

    
end

	break;
	
end
end

end
@label _out
begin
end
end
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

