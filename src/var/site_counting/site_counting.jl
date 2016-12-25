# Types and methods for counting mutations
# ========================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract SiteCase

abstract PairSiteCase <: SiteCase

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol.
"""
immutable Certain <: SiteCase end

"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
immutable Ambiguous <: SiteCase end

"""
An `Indel` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
immutable Indel <: SiteCase end

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
immutable Match <: PairSiteCase end

"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
immutable Mismatch <: PairSiteCase end

"""
A `Mismatch` site describes a site where two aligned nucleotides are definately
conserved. By definately conserved this means that the symbols of the site are
non-ambiguity symbols, and they are the same symbol.
"""
immutable Conserved <: PairSiteCase end

"""
A `Mutated` site describes a site where two aligned nucleotides are definately
mutated. By definately mutated this means that the symbols of the site are
non-ambiguity symbols, and they are not the same symbol.
"""
immutable Mutated <: PairSiteCase end

"""
A `Transition` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transition mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transition mutation: i.e. A<->G, or C<->T.
"""
immutable Transition <: PairSiteCase end

"""
A `Transversion` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transversion mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transversion mutation: i.e. A<->C,
A<->T, G<->T, G<->C.
"""
immutable Transversion <: PairSiteCase end

typealias FourBitAlphs Union{DNAAlphabet{4},RNAAlphabet{4}}
typealias TwoBitAlphs Union{DNAAlphabet{2},RNAAlphabet{2}}


# Includes for naieve site counting.
include("naieve/is_site.jl")
include("naieve/count_sites.jl")

#=

# Includes for performant site counting.
include("nibble_operations.jl")
include("nibble_counting.jl")
include("aligned_data_iterator.jl")



function count_sites(::Type{Ambiguous}, seq::BioSequence{FourBitAlphs})
    site_count = 0
    bindex = Seq.bitindex(seq, 1)
    firstoff = Seq.offset(bindex)

    # If the offset of the first element is not zero, then the first integer
    # needs to be shifted / masked.
    if firstoff != 0
        offsetint = seq.data[1] >> firstoff
        site_count += count_sites4(T, offsetint)
        for i in 2:endof(seq.data)
            site_count += count_sites4(Ambiguous, seq.data[i])
        end
    else
        for i in 1:endof(seq.data)
            site_count += count_sites4(Ambiguous, seq.data[i])
        end
    end
    return site_count
end


"""
    count_sites{T<:SiteCase,A}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})

Count the number of sites of tyle `T` between two BioSequences that are assumed
to be equal length and are therefore aligned.
"""
function count_sites{T<:SiteCase,A}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})

    if length(a) != length(b)
        error("BioSequences `a` and `b` are not of the same length.")
    end

    # First we check the offsets of A and B.
    # We check if the two sequences are subsequences or not, and then decide
    # what to do.
    oa = offset(bitindex(a, 1))
    ob = offset(bitindex(b, 1))

    ba = oa == 0
    bb = ob == 0

    site_count = 0

    if ba && bb
        for i in 1:endof(a.data)
            site_count += count_sites4(T, a.data[i], b.data[i])
        end
    elseif !ba && bb

    elseif


    return site_count
end








"""
    count_sites_loop{T<:SiteCase}(::Type{T}, a::Vector{UInt64}, b::Vector{UInt64})

Count the number of sites of a certain type in the binary data of two dna or
rna sequences.

This is an internal function, and it assumes:

- The encoding is the 4 bit encoding.
- That the binary data is aligned. i.e. that no bitshifting is nessecery since
  the two sequences are not subsequences or start in the middle of a vector.
"""
function count_sites_loop{T<:SiteCase}(::Type{T},
                                       a::Vector{UInt64},
                                       b::Vector{UInt64})
    n = 0
    @inbounds for i in eachindex(a)
        n += count_sites4(T, a[i], b[i])
    end
    return n
end

# Count the number of sites of a certain type in the binary data of two dna or rna sequences.
# Assume one of the sets of binary data is aligned, but the other is not, and so a ShiftedInts
# iterator is required.
function count_sites4{T<:SiteCase,A<:FourBitAlphs}(::Type{T}, a::ShiftedInts{A}, b::Vector{UInt64})
    if length(a) != length(b)
        error("`a` and `b` must be the same length.")
    end
    n = 0
    @inbounds for (i, j) in zip(a, b)
        n += count_sites4(T, i, j)
    end
    return n
end

function count_sites4{T<:SiteCase,A<:FourBitAlphs}(::Type{T}, a::Vector{UInt64}, b::ShiftedInts{A})
    return count_sites4(T, b, a)
end


@inline function get_shifted_integers(vec::Vector{UInt64}, state::BitIndex,
                                                     lastIndex::BitIndex,
                                                     lastMask::UInt64)
    firstIdx = index(state)

    # Determine if the state has reached the final integer containing sequence data.
    firstIsFinal = firstIdx == index(lastIndex)

    # Determine then if a second integer is needed from the data, and determine
    # if this second integer is the final integer containing sequence data.
    secondIdx = ifelse(firstIsFinal, index(state), index(state + 64))
    secondIsFinal = secondIdx == index(lastIndex)

    # Get the first and second integers from the vector.
    firstInt = vec[firstIdx]
    secondInt = vec[secondIdx]

    # If either of the two integers is the final integer containing sequence data,
    # we need to mask them to make sure the data really does end.
    firstInt &= ifelse(firstIsFinal, lastMask, 0xffffffffffffffff)
    secondInt &= ifelse(secondIsFinal, lastMask, 0xffffffffffffffff)

    firstInt = firstInt >> offset(state)

    firstInt |= ifelse(firstIsFinal, firstInt, secondInt << (64 - offset(state)))

    return firstInt
end

# Count the number of sites of a certain type in the binary data of two dna or rna sequences.
# Assume that one sequence needs shifting. In this case sequence 'a'.
function count_sites4{T<:SiteCase}(::Type{T},
                                   a::Vector{UInt64},
                                   b::Vector{UInt64},
                                   finalIndex::BitIndex,
                                   finalMask::UInt64,
                                   state::BitIndex,
                                   nIntegers::Int64)
    n = 0
    @inbounds for i in 1:nIntegers
        shiftedInt = get_shifted_integers(a, state, finalIndex, finalMask)
        n += count_sites4(T, shiftedInt, b[i])
        state += 64
    end
    return n
end

function count_sites4{T<:SiteCase,A<:FourBitAlphs}(::Type{T}, a::ShiftedInts{A}, b::ShiftedInts{A})
    @assert length(a) == length(b)
    n = 0
    @inbounds for (i, j) in zip(a, b)
        n += count_sites4(T, i, j)
    end
    return n
end

=#
