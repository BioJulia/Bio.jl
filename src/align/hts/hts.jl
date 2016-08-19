# HTS
# ===
#
# High-throughtput sequencing.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

import Bio: Ragel
import BGZFStreams: BGZFStream, VirtualOffset, virtualoffset
import BufferedStreams   # this is necessary though I don't know why.
import BufferedStreams: BufferedInputStream

include("chunk.jl")
include("sam/sam.jl")
include("bam/bam.jl")
include("tabix/tabix.jl")
