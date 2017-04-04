# site_types.jl
# =============
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
An abstract type which you can inherit from to build on the site counting
"""
abstract Site

# Methods required for both the naieve and bitparallel framework.
out_type{T<:Site}(::Type{T}) = Int
start_count{T<:Site}(::Type{T}) = 0

# Methods required for the naive framework.
update_count(acc::Int, up::Bool) = acc + up

# Methods required for the bitparallel framework.
update_count(acc::Int, up::Int) = acc + up
@inline correct_endspace{T<:Site}(::Type{T}) = false
@inline endspace_correction(nspace::Int, count::Int) = count - nspace

include("certain.jl")
include("ambiguous.jl")
include("gap.jl")
include("match.jl")
include("mismatch.jl")
