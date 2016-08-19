%%{
    machine samparser;

    action count_line {
        input.state.linenum += 1
    }

    action anchor {
        Ragel.@anchor!
    }

    action header {
        header = Ragel.@ascii_from_anchor!()
        tag = header[2:3]
        if !haskey(input.header, tag)
            input.header[tag] = []
        end
        if tag == "CO"
            push!(input.header[tag], header[5:end-1])
        else
            push!(input.header[tag], parse_keyvals(header[5:end-1]))
        end
    }

    action qname {
        Ragel.@copy_from_anchor!(output.name)
    }

    action flag {
        output.flag = Ragel.@int64_from_anchor!
    }

    action rname {
        Ragel.@copy_from_anchor!(output.refname)
    }

    action pos {
        output.pos = Ragel.@int64_from_anchor!
    }

    action mapq {
        output.mapq = Ragel.@int64_from_anchor!
    }

    action cigar {
        Ragel.@copy_from_anchor!(output.cigar)
    }

    action next_rname {
        Ragel.@copy_from_anchor!(output.next_refname)
    }

    action next_pos {
        output.next_pos = Ragel.@int64_from_anchor!
    }

    action tlen {
        output.tlen = Ragel.@int64_from_anchor!
    }

    action seq {
        seqstr = Ragel.@ascii_from_anchor!()
        resize!(output.seq, length(seqstr))
        Bio.Seq.encode_copy!(output.seq, seqstr)
    }

    action qual {
        qualstr = Ragel.@ascii_from_anchor!()
        resize!(output.qual, length(qualstr))
        for i in 1:endof(qualstr)
            output.qual[i] = UInt8(qualstr[i]) - 33
        end
        empty!(output.optional_fields)
    }

    action optfield {
        optfieldstr = Ragel.@ascii_from_anchor!()
        tag = optfieldstr[1:2]
        typ = optfieldstr[4]
        if typ == 'A'
            value = optfieldstr[6]
        elseif typ == 'Z'
            value = optfieldstr[6:end]
        elseif typ == 'H'
            value = parse_hexbytearray(optfieldstr[6:end])
        elseif typ == 'B'
            eltyp = auxtype[UInt8[optfieldstr[6]]]
            value = [parse(eltyp, x) for x in split(optfieldstr[8:end], ',')]
        else
            value = parse(auxtype[UInt8(typ)], optfieldstr[6:end])
        end
        output.optional_fields[tag] = value
    }

    action record {
        # need to anchor here for the next record
        Ragel.@anchor!
        Ragel.@yield ftargs
    }

    newline = '\n' >count_line;
    sign    = ('-' | '+')?;
    float   = sign? digit* '.'? digit+ ([eE] sign? digit+)?;
    tag     = alpha alnum;

    comment = (any - newline)*;
    keyval  = (tag ':' print+);

    header = ('@' ("CO\t" comment | tag ('\t' keyval)+) newline) >anchor %header;

    qname = (graph - '@')+ >anchor %qname;
    flag  = digit+ >anchor %flag;
    rname = graph+ >anchor %rname;
    pos   = digit+ >anchor %pos;
    mapq  = digit+ >anchor %mapq;
    cigar = ('*' | (digit+ [MIDNSHPX=])+) >anchor %cigar;
    next_rname = graph+ >anchor %next_rname;
    next_pos   = digit+ >anchor %next_pos;
    tlen  = ('-'? digit+) >anchor %tlen;
    seq   = ('*' | [A-Za-z=.]+) >anchor %seq;
    qual  = graph+ >anchor %qual;
    optfield = (tag ':' (
        ("A:" graph) |
        ("Z:" print+) |
        ("H:" [0-9A-F]+) |
        ("B:" (',' float)+) |
        ("f:" float) |
        ([cCsSiI] ':' sign? digit+))) >anchor %optfield;

    record = (
        qname '\t'
        flag  '\t'
        rname '\t'
        pos   '\t'
        mapq  '\t'
        cigar '\t'
        next_rname '\t'
        next_pos   '\t'
        tlen  '\t'
        seq   '\t'
        qual ('\t' optfield)*
        newline) %record;

    main := header* record*;
}%%

%% write data;

Ragel.@generate_read!_function(
    "samparser",
    SAMReader,
    SAMRecord,
    begin
        %% write exec;
    end)


# header only
%%{
    machine samheaderparser;

    action anchor {
        anchor = p + 1
    }

    action header {
        line = String(data[anchor:p])
        tag = line[2:3]
        if !haskey(header, tag)
            header[tag] = []
        end
        if tag == "CO"
            push!(header[tag], line[5:end-1])
        else
            push!(header[tag], parse_keyvals(line[5:end-1]))
        end
    }

    newline = '\n';
    tag     = alpha alnum;

    comment = (any - newline)*;
    keyval  = (tag ':' print+);

    header = ('@' ("CO\t" comment | tag ('\t' keyval)+) newline) >anchor %header;

    main := header*;
}%%

%% write data;

function parse_samheader(data)
    header = SAMHeader()

    p = 0
    pe = eof = endof(data)
    anchor = 0

    %% write init;
    %% write exec;

    @assert cs >= samheaderparser_first_final

    return header
end
