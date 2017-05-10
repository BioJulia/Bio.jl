# BBI
# ===
#
# The BigBed format is documented in
#   Kent, W. James, et al. "BigWig and BigBed: # enabling browsing of large
#   distributed datasets." Bioinformatics 26.17 (2010): 2204-2207.
#
# The low level details are documented in a series of tables in the supplement
# of that paper. The immutable types defined below are exactly the data layout
# described in those tables. They are labeled with the table they correspand to.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

include("types.jl")
include("reader.jl")
include("writer.jl")
include("intersection.jl")
