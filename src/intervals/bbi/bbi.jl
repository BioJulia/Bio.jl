# BBI - Shared Parts of bigWig and bigBed
# =======================================
#
# Detailed description of bigWig and bigBed file formats are found in the paper
# and its supplement:
# https://doi.org/10.1093/bioinformatics/btq351
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922891/bin/supp_btq351_bbiSuppFINAL.doc
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module BBI

import Bio: Mem
import Libz

include("header.jl")
include("summary.jl")
include("btree.jl")
include("rtree.jl")
include("section.jl")
include("zoom.jl")

function compress!(dst::Vector{UInt8}, src::Vector{UInt8})::UInt64
    len = Ref(Culong(length(dst)))
    code = ccall(
        (:compress, Libz.zlib),
        Cint,
        (Ptr{Void}, Ref{Culong}, Ptr{Void}, Culong),
        dst, len, src, sizeof(src))
    if code != Libz.Z_OK
        Libz.zerror(code)
    end
    return len[]
end

#=
function decompress()
end
=#

end
