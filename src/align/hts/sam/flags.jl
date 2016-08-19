# Flags
# =====
#
# Bitwise flags (or FLAG).

"the read is paired in sequencing, no matter whether it is mapped in a pair"
const SAM_FLAG_PAIRED        = UInt16(0x1)

"the read is mapped in a proper pair"
const SAM_FLAG_PROPER_PAIR   = UInt16(0x2)

"the read itself is unmapped; conflictive with SAM_FLAG_PROPER_PAIR"
const SAM_FLAG_UNMAP         = UInt16(0x4)

"the mate is unmapped"
const SAM_FLAG_MUNMAP        = UInt16(0x8)

"the read is mapped to the reverse strand"
const SAM_FLAG_REVERSE       = UInt16(0x10)

"the mate is mapped to the reverse strand"
const SAM_FLAG_MREVERSE      = UInt16(0x20)

"this is read1"
const SAM_FLAG_READ1         = UInt16(0x40)

"this is read2"
const SAM_FLAG_READ2         = UInt16(0x80)

"not primary alignment"
const SAM_FLAG_SECONDARY     = UInt16(0x100)

"QC failure"
const SAM_FLAG_QCFAIL        = UInt16(0x200)

"optical or PCR duplicate"
const SAM_FLAG_DUP           = UInt16(0x400)

"supplementary alignment"
const SAM_FLAG_SUPPLEMENTARY = UInt16(0x800)

