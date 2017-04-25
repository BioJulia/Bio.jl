# Extending BioSequence site counting to include mutations
# ========================================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract Mutation <: Bio.Seq.Site

Base.zero(::Type{Tuple{Int,Int}}) = zero(Int), zero(Int)

@inline function update_counter{T<:Union{Int,Bool}}(acc::Tuple{Int,Int}, up::Tuple{T,T})
    return (acc[1] + up[1]), (acc[2] + up[2])
end

@inline function issite{M<:Mutation}(::Type{M}, a::BioSequence, b::BioSequence, idx)
    k = ischange(M, a, b, idx)
    c = issite(Certain, a, b, idx)
    return k & c, c
end

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin
        counter_type{M<:Mutation}(::Type{M}, ::Type{$A{2}}, ::Type{$A{4}}) = Tuple{Int,Int}
        counter_type{M<:Mutation}(::Type{M}, ::Type{$A{4}}, ::Type{$A{2}}) = Tuple{Int,Int}
        counter_type{M<:Mutation}(::Type{M}, ::Type{$A{4}}, ::Type{$A{4}}) = Tuple{Int,Int}
        @inline function issite{M<:Mutation}(::Type{M}, a::BioSequence{$A{2}}, b::BioSequence{$A{2}}, idx)
            return ischange(M, a, b, idx)
        end
    end
end

include("conserved.jl")
include("mutated.jl")
include("transition.jl")
include("transversion.jl")
