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


const fastq_start  = convert(Int , 24)
const fastq_first_final  = convert(Int , 24)
const fastq_error  = convert(Int , 0)
const fastq_en_main  = convert(Int , 24)
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
        begin
cs = fastq_start;
	end
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
            begin
if p == pe 
	@goto _test_eof

end
@switch cs  begin
    @case 24
@goto st_case_24
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
@case 6
@goto st_case_6
@case 7
@goto st_case_7
@case 8
@goto st_case_8
@case 9
@goto st_case_9
@case 10
@goto st_case_10
@case 11
@goto st_case_11
@case 25
@goto st_case_25
@case 26
@goto st_case_26
@case 12
@goto st_case_12
@case 13
@goto st_case_13
@case 14
@goto st_case_14
@case 15
@goto st_case_15
@case 16
@goto st_case_16
@case 27
@goto st_case_27
@case 17
@goto st_case_17
@case 18
@goto st_case_18
@case 19
@goto st_case_19
@case 20
@goto st_case_20
@case 21
@goto st_case_21
@case 22
@goto st_case_22
@case 23
@goto st_case_23

end
@goto st_out
@label ctr0
begin
	input.state.linenum += 1
    
end
@goto st24
@label st24
p+= 1;
	if p == pe 
	@goto _test_eof24

end
@label st_case_24
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st24
end
@case 10
begin
@goto ctr0
end
@case 11
begin
@goto st24
end
@case 13
begin
@goto st1
end
@case 32
begin
@goto st24
end
@case 64
begin
ck  = convert(Int , 0)

if (length(input.qualbuf) == length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto st2
end

end
@goto st0
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
@label ctr52
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
        append!(input.namebuf, state.buffer, firstpos, p)
    
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
@label ctr12
begin
	input.state.linenum += 1
    
end
@goto st6
@label ctr5
begin
	firstpos = Ragel.@popmark!
        append!(input.namebuf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
@goto st6
@label ctr10
begin
	firstpos = Ragel.@popmark!
        append!(input.descbuf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
@goto st6
@label ctr41
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
@goto st6
@label st6
p+= 1;
	if p == pe 
	@goto _test_eof6

end
@label st_case_6
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr12
end
@case 13
begin
@goto st7
end
@case 43
begin
@goto st8
end

end
if ( data[1 + p ]) > 90 
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto ctr15
end

end
end

elseif ( ( data[1 + p ]) >= 65  )
	begin
@goto ctr15
end

end
begin
@goto st0
end
@label ctr6
begin
	firstpos = Ragel.@popmark!
        append!(input.namebuf, state.buffer, firstpos, p)
    
end
@goto st7
@label ctr11
begin
	firstpos = Ragel.@popmark!
        append!(input.descbuf, state.buffer, firstpos, p)
    
end
@goto st7
@label ctr42
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
@goto st7
@label st7
p+= 1;
	if p == pe 
	@goto _test_eof7

end
@label st_case_7
if ( data[1 + p ]) == 10 
	begin
@goto ctr12
end

end
begin
@goto st0
end
@label st8
p+= 1;
	if p == pe 
	@goto _test_eof8

end
@label st_case_8
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st8
end
@case 10
begin
@goto ctr17
end
@case 11
begin
@goto st8
end
@case 12
begin
@goto st0
end
@case 13
begin
@goto st13
end
@case 32
begin
@goto st8
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr16
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto ctr16
end

end
begin
@goto ctr16
end
@label ctr16
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
@goto ctr20
end
@case 10
begin
@goto ctr21
end
@case 11
begin
@goto ctr20
end
@case 12
begin
@goto st0
end
@case 13
begin
@goto ctr22
end
@case 32
begin
@goto ctr20
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto st9
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto st9
end

end
begin
@goto st9
end
@label ctr20
begin
	firstpos = Ragel.@popmark!
        append!(input.name2buf, state.buffer, firstpos, p)
    
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
@goto ctr24
end
@case 11
begin
@goto ctr24
end
@case 12
begin
@goto ctr23
end
@case 32
begin
@goto ctr24
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
@goto ctr23
end

end
end

else
	begin
@goto ctr23
end

end
begin
@goto ctr23
end
@label ctr23
begin
	Ragel.@pushmark!
    
end
@goto st11
@label st11
p+= 1;
	if p == pe 
	@goto _test_eof11

end
@label st_case_11
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr26
end
@case 13
begin
@goto ctr27
end

end
if ( data[1 + p ]) > 12 
	begin
if 14 <= ( data[1 + p ]) 
	begin
@goto st11
end

end
end

elseif ( ( data[1 + p ]) >= 11  )
	begin
@goto st11
end

end
begin
@goto st11
end
@label ctr17
begin
	input.state.linenum += 1
    
end
@goto st25
@label ctr21
begin
	firstpos = Ragel.@popmark!
        append!(input.name2buf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
@goto st25
@label ctr26
begin
	firstpos = Ragel.@popmark!
        append!(input.desc2buf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
@goto st25
@label ctr29
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
begin
	input.state.linenum += 1
    
end
@goto st25
@label ctr38
begin
	firstpos = Ragel.@popmark!
        append!(input.name2buf, state.buffer, firstpos, p)
    
end
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
begin
	input.state.linenum += 1
    
end
@goto st25
@label st25
p+= 1;
	if p == pe 
	@goto _test_eof25

end
@label st_case_25
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st26
end
@case 10
begin
@goto ctr17
end
@case 11
begin
@goto st26
end
@case 13
begin
@goto st13
end
@case 32
begin
@goto st26
end
@case 64
begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;

end
if (length(input.qualbuf) == length(input.seqbuf)
    )
	ck += 2;
	
end
if ck < 2 
	begin
if 1 <= ck 
	begin
@goto ctr51
end

end
end

elseif ( ck > 2  )
	begin
@goto ctr53
end

else
	begin
@goto ctr52
end

end
@goto st0
end

end
if ( data[1 + p ]) > 63 
	begin
if 65 <= ( data[1 + p ]) && ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr51
end

end
@goto st0
end

end
end

elseif ( ( data[1 + p ]) >= 33  )
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr51
end

end
@goto st0
end

end
begin
@goto st0
end
@label ctr28
begin
	input.state.linenum += 1
    
end
@goto st26
@label st26
p+= 1;
	if p == pe 
	@goto _test_eof26

end
@label st_case_26
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st26
end
@case 10
begin
@goto ctr28
end
@case 11
begin
@goto st26
end
@case 13
begin
@goto st12
end
@case 32
begin
@goto st26
end
@case 64
begin
ck  = convert(Int , 0)

if (length(input.qualbuf) == length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr52
end

end
@goto st0
end

end
begin
@goto st0
end
@label st12
p+= 1;
	if p == pe 
	@goto _test_eof12

end
@label st_case_12
if ( data[1 + p ]) == 10 
	begin
@goto ctr28
end

end
begin
@goto st0
end
@label ctr22
begin
	firstpos = Ragel.@popmark!
        append!(input.name2buf, state.buffer, firstpos, p)
    
end
@goto st13
@label ctr27
begin
	firstpos = Ragel.@popmark!
        append!(input.desc2buf, state.buffer, firstpos, p)
    
end
@goto st13
@label ctr30
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
@goto st13
@label ctr39
begin
	firstpos = Ragel.@popmark!
        append!(input.name2buf, state.buffer, firstpos, p)
    
end
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
@goto st13
@label st13
p+= 1;
	if p == pe 
	@goto _test_eof13

end
@label st_case_13
if ( data[1 + p ]) == 10 
	begin
@goto ctr17
end

end
begin
@goto st0
end
@label ctr31
begin
	input.qualcount += 1
    
end
@goto st14
@label ctr51
begin
	Ragel.@pushmark!
    
end
begin
	input.qualcount += 1
    
end
@goto st14
@label st14
p+= 1;
	if p == pe 
	@goto _test_eof14

end
@label st_case_14
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr29
end
@case 13
begin
@goto ctr30
end

end
if 33 <= ( data[1 + p ]) && ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr31
end

end
@goto st0
end

end
begin
@goto st0
end
@label ctr53
begin
	Ragel.@pushmark!
    
end
begin
	input.qualcount += 1
    
end
begin
	yield = true;
        begin
	p+= 1; cs = 15; @goto _out

end

    
end
@goto st15
@label st15
p+= 1;
	if p == pe 
	@goto _test_eof15

end
@label st_case_15
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr29
end
@case 13
begin
@goto ctr30
end
@case 32
begin
@goto st0
end
@case 127
begin
@goto ctr2
end

end
if ( data[1 + p ]) < 14 
	begin
if 9 <= ( data[1 + p ]) && ( data[1 + p ]) <= 12 
	begin
@goto st0
end

end
end

elseif ( ( data[1 + p ]) > 31  )
	begin
if 33 <= ( data[1 + p ]) && ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr32
end

else
	begin
@goto ctr2
end

end
@goto st0
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
@label ctr35
begin
	input.qualcount += 1
    
end
@goto st16
@label ctr32
begin
	Ragel.@pushmark!
    
end
begin
	input.qualcount += 1
    
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
@goto ctr4
end
@case 10
begin
@goto ctr33
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
@goto ctr34
end
@case 32
begin
@goto ctr4
end
@case 127
begin
@goto st3
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) && ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr35
end

else
	begin
@goto st3
end

end
@goto st0
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
@label ctr36
begin
	input.state.linenum += 1
    
end
@goto st27
@label ctr33
begin
	firstpos = Ragel.@popmark!
        append!(input.namebuf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
@goto st27
@label ctr44
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
begin
	input.state.linenum += 1
    
end
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
@goto st27
@label st27
p+= 1;
	if p == pe 
	@goto _test_eof27

end
@label st_case_27
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st26
end
@case 10
begin
@goto ctr36
end
@case 11
begin
@goto st26
end
@case 13
begin
@goto st17
end
@case 32
begin
@goto st26
end
@case 43
begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr56
end

else
	begin
@goto st8
end

end
@goto st0
end
@case 64
begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;

end
if (length(input.qualbuf) == length(input.seqbuf)
    )
	ck += 2;
	
end
if ck < 2 
	begin
if 1 <= ck 
	begin
@goto ctr51
end

end
end

elseif ( ck > 2  )
	begin
@goto ctr53
end

else
	begin
@goto ctr52
end

end
@goto st0
end

end
if ( data[1 + p ]) < 65 
	begin
if ( data[1 + p ]) > 42 
	begin
if 44 <= ( data[1 + p ]) && ( data[1 + p ]) <= 63 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr51
end

end
@goto st0
end

end
end

elseif ( ( data[1 + p ]) >= 33  )
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr51
end

end
@goto st0
end

end
end

elseif ( ( data[1 + p ]) > 90  )
	begin
if ( data[1 + p ]) < 97 
	begin
begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr51
end

end
@goto st0
end
end

elseif ( ( data[1 + p ]) > 122  )
	begin
if ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr51
end

end
@goto st0
end

end
end

else
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr57
end

else
	begin
@goto ctr15
end

end
@goto st0
end

end
end

else
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr57
end

else
	begin
@goto ctr15
end

end
@goto st0
end

end
begin
@goto st0
end
@label ctr34
begin
	firstpos = Ragel.@popmark!
        append!(input.namebuf, state.buffer, firstpos, p)
    
end
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
@goto st17
@label ctr45
begin
	firstpos = Ragel.@popmark!
        append!(input.seqbuf, state.buffer, firstpos, p)
    
end
begin
	firstpos = Ragel.@popmark!
        append!(input.qualbuf, state.buffer, firstpos, p)
        input.qualcount = 0
    
end
@goto st17
@label st17
p+= 1;
	if p == pe 
	@goto _test_eof17

end
@label st_case_17
if ( data[1 + p ]) == 10 
	begin
@goto ctr36
end

end
begin
@goto st0
end
@label ctr56
begin
	Ragel.@pushmark!
    
end
begin
	input.qualcount += 1
    
end
@goto st18
@label st18
p+= 1;
	if p == pe 
	@goto _test_eof18

end
@label st_case_18
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto st8
end
@case 10
begin
@goto ctr29
end
@case 11
begin
@goto st8
end
@case 12
begin
@goto st0
end
@case 13
begin
@goto ctr30
end
@case 32
begin
@goto st8
end
@case 127
begin
@goto ctr16
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) && ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr37
end

else
	begin
@goto ctr16
end

end
@goto st0
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto ctr16
end

end
begin
@goto ctr16
end
@label ctr40
begin
	input.qualcount += 1
    
end
@goto st19
@label ctr37
begin
	Ragel.@pushmark!
    
end
begin
	input.qualcount += 1
    
end
@goto st19
@label st19
p+= 1;
	if p == pe 
	@goto _test_eof19

end
@label st_case_19
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr20
end
@case 10
begin
@goto ctr38
end
@case 11
begin
@goto ctr20
end
@case 12
begin
@goto st0
end
@case 13
begin
@goto ctr39
end
@case 32
begin
@goto ctr20
end
@case 127
begin
@goto st9
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) && ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr40
end

else
	begin
@goto st9
end

end
@goto st0
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto st9
end

end
begin
@goto st9
end
@label ctr15
begin
	Ragel.@pushmark!
    
end
@goto st20
@label st20
p+= 1;
	if p == pe 
	@goto _test_eof20

end
@label st_case_20
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr41
end
@case 13
begin
@goto ctr42
end

end
if ( data[1 + p ]) > 90 
	begin
if 97 <= ( data[1 + p ]) && ( data[1 + p ]) <= 122 
	begin
@goto st20
end

end
end

elseif ( ( data[1 + p ]) >= 65  )
	begin
@goto st20
end

end
begin
@goto st0
end
@label ctr46
begin
	input.qualcount += 1
    
end
@goto st21
@label ctr57
begin
	Ragel.@pushmark!
    
end
begin
	Ragel.@pushmark!
    
end
begin
	input.qualcount += 1
    
end
@goto st21
@label st21
p+= 1;
	if p == pe 
	@goto _test_eof21

end
@label st_case_21
@switch ( data[1 + p ])  begin
    @case 10
begin
@goto ctr44
end
@case 13
begin
@goto ctr45
end

end
if ( data[1 + p ]) < 91 
	begin
if ( data[1 + p ]) > 64 
	begin
begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr46
end

else
	begin
@goto st20
end

end
@goto st0
end
end

elseif ( ( data[1 + p ]) >= 33  )
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr31
end

end
@goto st0
end

end
end

elseif ( ( data[1 + p ]) > 96  )
	begin
if ( data[1 + p ]) > 122 
	begin
if ( data[1 + p ]) <= 126 
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr31
end

end
@goto st0
end

end
end

else
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if ck > 0 
	begin
@goto ctr46
end

else
	begin
@goto st20
end

end
@goto st0
end

end
end

else
	begin
ck  = convert(Int , 0)

if (length(input.qualbuf) + input.qualcount < length(input.seqbuf)
    )
	ck += 1;
	
end
if 1 <= ck 
	begin
@goto ctr31
end

end
@goto st0
end

end
begin
@goto st0
end
@label ctr24
begin
	Ragel.@pushmark!
    
end
@goto st22
@label st22
p+= 1;
	if p == pe 
	@goto _test_eof22

end
@label st_case_22
@switch ( data[1 + p ])  begin
    @case 9
begin
@goto ctr24
end
@case 10
begin
@goto ctr26
end
@case 11
begin
@goto ctr24
end
@case 12
begin
@goto ctr23
end
@case 13
begin
@goto ctr27
end
@case 32
begin
@goto ctr24
end

end
if ( data[1 + p ]) > 31 
	begin
if 33 <= ( data[1 + p ]) 
	begin
@goto ctr23
end

end
end

elseif ( ( data[1 + p ]) >= 14  )
	begin
@goto ctr23
end

end
begin
@goto ctr23
end
@label ctr8
begin
	Ragel.@pushmark!
    
end
@goto st23
@label st23
p+= 1;
	if p == pe 
	@goto _test_eof23

end
@label st_case_23
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
@label _test_eof24
cs = 24; @goto _test_eof
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
@label _test_eof6
cs = 6; @goto _test_eof
@label _test_eof7
cs = 7; @goto _test_eof
@label _test_eof8
cs = 8; @goto _test_eof
@label _test_eof9
cs = 9; @goto _test_eof
@label _test_eof10
cs = 10; @goto _test_eof
@label _test_eof11
cs = 11; @goto _test_eof
@label _test_eof25
cs = 25; @goto _test_eof
@label _test_eof26
cs = 26; @goto _test_eof
@label _test_eof12
cs = 12; @goto _test_eof
@label _test_eof13
cs = 13; @goto _test_eof
@label _test_eof14
cs = 14; @goto _test_eof
@label _test_eof15
cs = 15; @goto _test_eof
@label _test_eof16
cs = 16; @goto _test_eof
@label _test_eof27
cs = 27; @goto _test_eof
@label _test_eof17
cs = 17; @goto _test_eof
@label _test_eof18
cs = 18; @goto _test_eof
@label _test_eof19
cs = 19; @goto _test_eof
@label _test_eof20
cs = 20; @goto _test_eof
@label _test_eof21
cs = 21; @goto _test_eof
@label _test_eof22
cs = 22; @goto _test_eof
@label _test_eof23
cs = 23; @goto _test_eof
@label _test_eof
begin
end
if p == eof 
	begin
@switch cs  begin
    @case 25
@case 26
@case 27
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


