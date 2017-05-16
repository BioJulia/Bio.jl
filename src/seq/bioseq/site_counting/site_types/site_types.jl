# site_types.jl
# =============
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md




# For a given site type you must define the following functions for the
# generated function and method dispatch to work on:

# Methods required for both the naive and bitparallel algorithms.
# ---------------------------------------------------------------

# What will the output type of count(YOUR_SITE_TYPE, seqa, seqb), be?
# Note that the type specified here, should have a Base.zero or start_counter
# method defined for it.
