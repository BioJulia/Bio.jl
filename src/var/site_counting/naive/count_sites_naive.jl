# Methods for site counting in sequences using a naive method
# ============================================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Methods for a single sequence.

"""
    count_sites_naive{T<:Site,A<:FourBitAlphs}(::Type{T}, a::BioSequence{A})
"""
function count_sites_naive{T<:Site,A<:FourBitAlphs}(::Type{T}, a::BioSequence{A})
    k = 0
    @inbounds for idx in eachindex(a)
        k += issite(T, a[idx])
    end
    return k
end

count_sites_naive{T<:Union{Gap,Ambiguous},A<:TwoBitAlphs}(::Type{T}, a::BioSequence{A}) = 0
count_sites_naive{A<:TwoBitAlphs}(::Type{Certain}, a::BioSequence{A}) = length(a)

"""
    count_sites_naive{T<:Site,A<:DNA_OR_RNA}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})
This method counts the number of sites between a pair of aligned sequences of
type `T`.
"""
function count_sites_naive{T<:Site,A<:Alphabet,B<:Alphabet}(::Type{T},
                                                       a::BioSequence{A},
                                                       b::BioSequence{B})
    k = 0
    @inbounds for idx in 1:min(endof(a), endof(b))
        k += issite(T, a[idx], b[idx])
    end
    return k
end

"""
    count_sites_naive{T<:Mutation,A<:DNA_OR_RNA}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})
This method counts the number of sites between a pair of aligned sequences of
SiteCase `T`.
Since the types of SiteCase counted by this function may not always be
unambiguously determined, this method returns a tuple of the count as well as
the number of sites which could not be unambiguously determined.
This second count is important for some downstream purposes, for example
evolutionary/genetic distance computations in which pairwise deletion of
ambiguous sites is nessecery.
"""
function count_sites_naive{T<:Mutation,A<:Alphabet,B<:Alphabet}(::Type{T},
                                                             a::BioSequence{A},
                                                             b::BioSequence{B})
    k = 0
    j = 0
    @inbounds for idx in 1:min(endof(a), endof(b))
        an = a[idx]
        bn = b[idx]
        k += issite(T, an, bn)
        j += !issite(Certain, an, bn)
    end
    return k, j
end

count_sites_naive{T<:Union{Gap,Ambiguous},A<:TwoBitAlphs}(::Type{T}, a::BioSequence{A}, b::BioSequence{A}) = 0
count_sites_naive{A<:TwoBitAlphs}(::Type{Certain}, a::BioSequence{A}, b::BioSequence{A}) = min(length(a), length(b))
