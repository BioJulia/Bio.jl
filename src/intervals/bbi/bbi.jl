# BBI - Shared Parts of bigWig and bigBed
# =======================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module BBI

import Bio: Mem

include("header.jl")
include("summary.jl")
include("btree.jl")
include("rtree.jl")
include("section.jl")
include("zoom.jl")

end
