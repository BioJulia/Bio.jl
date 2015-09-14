
using Bio: FileFormat, AbstractParser
using Bio.StringFields
using Bio.Ragel
using BufferedStreams
using Switch
import Base: open, read!, eltype, copy

immutable BED3 <: FileFormat end

type BED3Entry
    seqname::StringField
    start::Int64
    stop::Int64
end


function BED3Entry()
    return BED3Entry(StringField(), 1, 0)
end


function copy(entry::BED3Entry)
    return BED3Entry(copy(entry.seqname), entry.start, entry.stop)
end


%%{
    machine bed3;

    action anchor {
        Ragel.@anchor!
    }

    action finish_match {
        # calling anchor is necessary here because fbreak will skip the entry
        # action for the next state (seqname, in this case)
        Ragel.@anchor!
        yield = true;
        fbreak;
    }

    action seqname {
        Ragel.@copy_from_anchor!(output.seqname)
    }

    action start {
        # convert from 0-based to 1-based
        output.start = 1 + Ragel.@int64_from_anchor!
    }

    action stop {
        output.stop = Ragel.@int64_from_anchor!
    }

    seqname = [ -~]*   >anchor   %seqname;
    start   = digit+   >anchor   %start;
    stop    = digit+   >anchor   %stop;
    bed3_entry = seqname '\t' start '\t' stop '\n';
    main := (bed3_entry %finish_match)*;
}%%


%% write data;


type BED3Parser <: AbstractParser
    state::Ragel.State
end


function eltype(::Type{BED3Parser})
    return BED3Entry
end


function open(input::BufferedInputStream, ::Type{BED3})
    %% write init;

    return BED3Parser(Ragel.State(cs, input))
end


Ragel.@generate_read_fuction("bed3", BED3Parser, BED3Entry,
    begin
        %% write exec;
    end)


