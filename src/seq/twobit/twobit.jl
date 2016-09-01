# 2bit format
# ===========
#
# Reader and writer of the 2bit file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
The 2bit file format stores multiple DNA sequences (up to 4 Gbp total) as a
compact randomly-accessible format.
See https://genome.ucsc.edu/FAQ/FAQformat.html#format7 for the details.
"""
immutable TwoBit <: Bio.IO.FileFormat end

include("reader.jl")
include("writer.jl")
