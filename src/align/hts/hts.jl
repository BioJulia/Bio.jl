# HTS
# ===
#
# High-throughtput sequencing.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

include("sam/sam.jl")
include("bam/bam.jl")

function Base.show(io::IO, rec::Union{SAMRecord,BAMRecord})
    name = seqname(rec)
    if isempty(name)
        name = "(empty name)"
    end
    if get(io, :compact, false)
        if ismapped(rec)
            range = string(refname(rec), ':', leftposition(rec), '-', rightposition(rec))
            print(io, name, "\t$(range)\t", cigar(rec))
        else
            print(io, name, "\tunmapped")
        end
    else
        println(io, summary(rec), ':')
        println(io, "  sequence name: ", name)
        println(io, "  reference name: ", refname(rec))
        println(io, "  leftmost position: ", leftposition(rec))
        println(io, "  next reference name: ", nextrefname(rec))
        println(io, "  next leftmost position: ", nextleftposition(rec))
        println(io, "  mapping quality: ", mappingquality(rec))
        println(io, "  flag: ", flag(rec))
        println(io, "  template length: ", templatelength(rec))
        println(io, "  CIGAR string: ", cigar(rec))
        println(io, "  sequence: ", String(sequence(rec)))
        println(io, "  base qualities: ", qualities(rec))
          print(io, "  optional fields: ", optional_fields(rec))
    end
end
