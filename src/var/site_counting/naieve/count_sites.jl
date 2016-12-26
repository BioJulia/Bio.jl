# Methods for site counting in sequences using a naieve method
# ============================================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Methods for a single sequence.

"""
    count_sites_naieve{T<:SiteCase,A<:FourBitAlphs}(::Type{T}, a::BioSequence{A})


"""
function count_sites_naieve{T<:SiteCase,A<:FourBitAlphs}(::Type{T}, a::BioSequence{A})
    k = 0
    @inbounds for idx in eachindex(a)
        k += ifelse(issite(T, a[idx]), 1, 0)
    end
    return k
end

count_sites_naieve{T<:Union{Indel,Ambiguous},A<:TwoBitAlphs}(::Type{T}, a::BioSequence{A}) = 0
count_sites_naieve{A<:TwoBitAlphs}(::Type{Certain}, a::BioSequence{A}) = length(a)

"""
    count_sites_naieve{T<:SiteCase,A<:Alphabet,B<:Alphabet}(::Type{T}, a::BioSequence{A}, b::BioSequence{B})

This method counts the number of sites between a pair of aligned sequences of
type `T`.
"""
function count_sites_naieve{T<:SiteCase,A<:Alphabet,B<:Alphabet}(::Type{T},
                                                                 a::BioSequence{A},
                                                                 b::BioSequence{B})
    k = 0
    @inbounds for idx in 1:min(endof(a), endof(b))
        k += ifelse(issite(T, a[idx], b[idx]), 1, 0)
    end
    return k
end

"""
    count_sites_naieve{T<:PairSiteCase,A<:Alphabet,B<:Alphabet}(::Type{T}, a::BioSequence{A}, b::BioSequence{B})

This method counts the number of sites between a pair of aligned sequences of
SiteCase `T`.

Since the types of SiteCase counted by this function may not always be
unambiguously determined, this method returns a tuple of the count as well as
the number of sites which could not be unambiguously determined.

This second count is important for some downstream purposes, for example
evolutionary/genetic distance computations in which pairwise deletion of
ambiguous sites is nessecery.
"""
function count_sites_naieve{T<:PairSiteCase,A<:Alphabet,B<:Alphabet}(::Type{T},
                                                                     a::BioSequence{A},
                                                                     b::BioSequence{B})
    k = 0
    j = 0
    @inbounds for idx in 1:min(endof(a), endof(b))
        an = a[idx]
        bn = b[idx]
        k += issite(T, an, bn)
        j += (issite(Ambiguous, an, bn) | issite(Indel, an, bn))
    end
    return k, j
end

count_sites_naieve{T<:Union{Indel,Ambiguous},A<:TwoBitAlphs}(::Type{T}, a::BioSequence{A}, b::BioSequence{A}) = 0
count_sites_naieve{A<:TwoBitAlphs}(::Type{Certain}, a::BioSequence{A}, b::BioSequence{A}) = min(length(a), length(b))
