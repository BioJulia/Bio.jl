# Types and methods for counting mutations
# ========================================
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract SiteCase
immutable Conserved <: SiteCase end
immutable Mutated <: SiteCase end
immutable Transition <: SiteCase end
immutable Transversion <: SiteCase end
immutable Gap <: SiteCase end
immutable Ambiguous <: SiteCase end
immutable Pairdel <: SiteCase end

include("bitwise_ops.jl")
include("aligned_data_iterator.jl")

typealias FourBitAlphs Union{DNAAlphabet{4},RNAAlphabet{4}}

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


# Count the number of sites of a certain type in the binary data of two dna or
# rna sequences. Assume the encoding is the 4 bit encoding. Also assume that the
# binary data is aligned.
function count_sites4{T<:SiteCase}(::Type{T}, a::Vector{UInt64}, b::Vector{UInt64})
    if length(a) != length(b)
        error("'a' and 'b' must be the same length.")
    end
    n = 0
    @inbounds for i in 1:endof(a)
        n += count_sites4(T, a[i], b[i])
    end
    return n
end

# Count the number of sites of a certain type in the binary data of two dna or rna sequences.
# Assume one of the sets of binary data is aligned, but the other is not, and so a ShiftedInts
# iterator is required.
function count_sites4{T<:SiteCase,A<:FourBitAlphs}(::Type{T}, a::ShiftedInts{A}, b::Vector{UInt64})
    if length(a) != length(b)
        error("'a' and 'b' must be the same length.")
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

#=
function count_sites{T<:SiteCase,A}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})
    @assert length(a) == length(b)

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
=#
