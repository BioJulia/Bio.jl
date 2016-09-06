# Index
# =====
#
# Index types for genomic intervals.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

import BGZFStreams: BGZFStream, VirtualOffset

include("chunk.jl")
include("bgzfindex.jl")
include("tabix.jl")
