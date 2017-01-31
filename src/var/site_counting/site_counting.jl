# Types and methods for counting mutations
# ========================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
The type `SiteCase` is a parametric abstract type that is used to organise the
different types of site that are available for counting in the Var site counting
framework.
Each type in the Var site counting framework inhertis from either
`SiteCase{true}` or `SiteCase{false}`.
If the type inhertis from `SiteCase{true}`, this means that it is not always
possible to unambiguously determine if a site is or isn't of the given type.
As an example, consider a site where one sequence has an an `A` nucleotide and
the other sequence has an `M` ambiguity code: the site could be a site of type
`Mutation`, but also could be a site of type `Conserved`, depending on what `M`
actually is in reality.
As a result, counting methods for such types that inherit from `SiteCase{true}`
also count and return the number of sites which could not be unambiguously
determined.
This second count is important for some downstream purposes, for example
evolutionary/genetic distance computations in which pairwise deletion of
ambiguous sites is necessary.
"""
abstract SiteCase{p}

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol.
"""
immutable Certain <: SiteCase{false} end

"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
immutable Ambiguous <: SiteCase{false} end

"""
An `Indel` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
immutable Indel <: SiteCase{false} end

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
immutable Match <: SiteCase{false} end

"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
immutable Mismatch <: SiteCase{false} end

"""
A `Mismatch` site describes a site where two aligned nucleotides are definately
conserved. By definately conserved this means that the symbols of the site are
non-ambiguity symbols, and they are the same symbol.
"""
immutable Conserved <: SiteCase{true} end

"""
A `Mutated` site describes a site where two aligned nucleotides are definately
mutated. By definately mutated this means that the symbols of the site are
non-ambiguity symbols, and they are not the same symbol.
"""
immutable Mutated <: SiteCase{true} end

"""
A `Transition` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transition mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transition mutation: i.e. A<->G, or C<->T.
"""
immutable Transition <: SiteCase{true} end

"""
A `Transversion` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transversion mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transversion mutation: i.e. A<->C,
A<->T, G<->T, G<->C.
"""
immutable Transversion <: SiteCase{true} end

typealias FourBitAlphs Union{DNAAlphabet{4},RNAAlphabet{4}}
typealias TwoBitAlphs Union{DNAAlphabet{2},RNAAlphabet{2}}

# Includes for naieve site counting.
include("naive/is_site.jl")
include("naive/count_sites_naive.jl")

"""
    count_sites{T<:SiteCase}(::Type{T}, a::BioSequence, b::BioSequence)

Mutations are counted using the `count_sites` method.

If the concrete site case type inherits from SiteCase{false}, then the result will
be a simple integer value of the number of counts of site.

If the concrete site case type inherits from SiteCase{true}, then the result
will be a tuple. As a result, counting methods for such types also count and
return the number of sites which could not be unambiguously determined.
This second count is important for some downstream purposes, for example
evolutionary/genetic distance computations in which pairwise deletion of
ambiguous sites is necessary.
"""
function count_sites{T<:SiteCase}(::Type{T}, a::BioSequence, b::BioSequence)
    return count_sites_naive(T, a, b)
end

function count_sites{T<:SiteCase}(::Type{T}, seq::BioSequence)
    return count_sites_naive(T, seq)
end
