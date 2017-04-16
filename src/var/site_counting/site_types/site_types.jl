# Extending BioSequence site counting to include mutations
# ========================================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

Base.zero(::Type{Tuple{Int,Int}}) = zero(Int), zero(Int)

abstract Mutation <: Site

counter_type{M<:Mutation,A<:Alphabet,B<:Alphabet}(::Type{M}, ::Type{A}, ::Type{B}) = Tuple{Int,Int}

update_counter(kacc::Int, cacc::Int, k::Bool, c::Bool) = kacc + k, cacc + c

@inline function issite{T<:NucleicAcid,M<:Mutation}(::Type{M}, a::T, b::T)
    k = ischange(M, a, b)
    c = issite(Certain, a, b)
    return k & c, c
end

@inline function count_nibbles{}()

end

include("conserved.jl")
include("mutated.jl")
include("transition.jl")
include("transversion.jl")
