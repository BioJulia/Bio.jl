# HTS
# ===
#
# High-throughtput sequencing.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

import Bio: Ragel
import Bio.Intervals: Chunk, BGZFIndex, BinIndex, LinearIndex, PseudoBin, LinearWindowSize
import BGZFStreams: BGZFStream, VirtualOffset, virtualoffset
import BufferedStreams   # this is necessary though I don't know why.
import BufferedStreams: BufferedInputStream

include("sam/sam.jl")
include("bam/bam.jl")
