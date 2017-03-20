# Flags
# =====
#
# Bitwise flags (or FLAG).

import Base: enumerate

baremodule SAMFlags

    immutable SAMFlag
        flag::UInt16
        description::String
    end
    
    const PAIRED        = SAMFlag(UInt16(0x1),
    "the read is paired in sequencing, no matter whether it is mapped in a pair")
    
    const PROPER_PAIR   = SAMFlag(UInt16(0x2),
    "the read is mapped in a proper pair")

    const UNMAP         = SAMFlag(UInt16(0x4),
    "the read itself is unmapped; conflictive with PROPER_PAIR")

    const MUNMAP        = SAMFlag(UInt16(0x8),
    "the mate is unmapped")

    const REVERSE       = SAMFlag(UInt16(0x10),
    "the read is mapped to the reverse strand")

    const MREVERSE      = SAMFlag(UInt16(0x20),
    "the mate is mapped to the reverse strand")

    const READ1         = SAMFlag(UInt16(0x40),
    "this is read1 (first in pair)")

    const READ2         = SAMFlag(UInt16(0x80),
    "this is read2 (second in pair)")

    const SECONDARY     = SAMFlag(UInt16(0x100),
    "not primary alignment")

    const QCFAIL        = SAMFlag(UInt16(0x200),
    "QC failure")

    const DUP           = SAMFlag(UInt16(0x400),
    "optical or PCR duplicate")

    const SUPPLEMENTARY = SAMFlag(UInt16(0x800),
    "supplementary alignment")
end

flag(f::SAMFlags.SAMFlag) = f.flag

function enumerate(::Type{SAMFlags.SAMFlag})
    n = names(SAMFlags,true)
    v = map(x->getfield(SAMFlags,x),n)
    
    out = SAMFlags.SAMFlag[]
    for x in v 
        typeof(x) == SAMFlags.SAMFlag && push!(out,x)
    end
    out
end

isflag(rec::Union{Bio.Align.SAMRecord,Bio.Align.BAMRecord},f::SAMFlags.SAMFlag) = (flag(rec) & flag(f)) == flag(f)

function decodeflag(rec::Union{Bio.Align.SAMRecord,Bio.Align.BAMRecord})
    for f in enumerate(SAMFlags.SAMFlag)
        if isflag(rec,f)
            println(f.description)
        end
    end
end
