# SAM Flags
# =========
#
# Bitwise flags (or FLAG).

"the read is paired in sequencing, no matter whether it is mapped in a pair"
const FLAG_PAIRED        = UInt16(0x1)

"the read is mapped in a proper pair"
const FLAG_PROPER_PAIR   = UInt16(0x2)

"the read itself is unmapped; conflictive with SAM.FLAG_PROPER_PAIR"
const FLAG_UNMAP         = UInt16(0x4)

"the mate is unmapped"
const FLAG_MUNMAP        = UInt16(0x8)

"the read is mapped to the reverse strand"
const FLAG_REVERSE       = UInt16(0x10)

"the mate is mapped to the reverse strand"
const FLAG_MREVERSE      = UInt16(0x20)

"this is read1"
const FLAG_READ1         = UInt16(0x40)

"this is read2"
const FLAG_READ2         = UInt16(0x80)

"not primary alignment"
const FLAG_SECONDARY     = UInt16(0x100)

"QC failure"
const FLAG_QCFAIL        = UInt16(0x200)

"optical or PCR duplicate"
const FLAG_DUP           = UInt16(0x400)

"supplementary alignment"
const FLAG_SUPPLEMENTARY = UInt16(0x800)
