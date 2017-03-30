# site_counting.jl
# ================
#
# A site counting framework for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

include("site_types.jl")
include("count_sites_naive.jl")

function count{T<:Site}(::Type{T}, a::BioSequence)
    return count_sites_naieve(T, a)
end

function count{T<:Site}(::Type{T}, a::BioSequence, b::BioSequence)
    return count_sites_naive(T, a, b)
end
