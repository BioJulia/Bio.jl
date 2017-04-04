# site_types.jl
# =============
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
An abstract type which you can inherit from to build on the site counting
"""
abstract Site
out_type{T<:Site}(::Type{T}) = Int
start_count{T<:Site}(::Type{T}) = 0
update_count(acc::Int, up::Bool) = acc + up


# Define Certain sites for the site-counting framework.
# -----------------------------------------------------

## The type.
"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol or a gap.
"""
immutable Certain <: Site end

## Methods for the naive framework.
"Test whether a nucleic acid in a sequence is certain."
issite{T<:NucleicAcid}(::Type{Certain}, a::T) = iscertain(a)

"Test whether a nucleotide site of two aligned sequences has two certain NucleicAcids."
@inline function issite{T<:NucleicAcid}(::Type{Certain}, a::T, b::T)
    return issite(Certain, a) & issite(Certain, b)
end


# Define Ambiguous sites for the site-counting framework.
# -------------------------------------------------------

## The type.
"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
immutable Ambiguous <: Site end

## Methods for the naive framework.
"Test whether a nucleic acid in a sequence is ambiguous."
issite{T<:NucleicAcid}(::Type{Ambiguous}, a::T) = isambiguous(a)

"Test whether a nucleotide site of two aligned sequences has ambiguities."
@inline function issite{T<:NucleicAcid}(::Type{Ambiguous}, a::T, b::T)
    return issite(Ambiguous, a) | issite(Ambiguous, b)
end


# Define Gaps for the site-counting framework.
# --------------------------------------------

## The type.
"""
An `Gap` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
immutable Gap <: Site end

## Methods for naive framework.
"Test whether a nucleic acid in a sequence is a gap character."
issite{T<:NucleicAcid}(::Type{Gap}, a::T) = isgap(a)

"Test whether a nucleotide site of two aligned sequences has gap characters."
@inline function issite{T<:NucleicAcid}(::Type{Gap}, a::T, b::T)
    return issite(Gap, a) | issite(Gap, b)
end


# Define Matches for the site-counting framework.
# -----------------------------------------------

## The type.
"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
immutable Match <: Site end

## Methods for the naive framework.
"Test whether a nucleotide site of two aligned sequences, constitutes a match."
issite{T<:NucleicAcid}(::Type{Match}, a::T, b::T) = a == b


# Define mismatches for site-counting framework.
# ----------------------------------------------

## The type.
"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
immutable Mismatch <: Site end

## Methods for the naive framework.
"Test whether a nucleotide site of two aligned sequences, constitutes a mismatch."
issite{T<:NucleicAcid}(::Type{Mismatch}, a::T, b::T) = a != b
