

module BED3ParserImpl

%%{
    machine bed3;

    action mark {
        Ragel.@mark!
    }

    action seqname {
        state.seqname = Ragel.@asciistring_from_mark!
    }

    action start {
        state.start = Ragel.@int64_from_mark!
    }

    action stop {
        state.stop = Ragel.@int64_from_mark!
    }

    seqname = [ -~]*   >mark   %seqname;
    start   = digit+   >mark   %start;
    stop    = digit+   >mark   %stop;
    bed3_entry = seqname '\t' start '\t' stop '\n';
    main := bed3_entry*;
}%%



end # module BED3ParserImpl