
%%{
    machine bed3;

    seqname = [ -~]*;
    first   = digit+;
    last    = digit+;
    bed3_entry = seqname '\t' first '\t' last '\n';
    main := bed3_entry*;
}%%

