# Iterating over binary sequence data in an aligned manner
# ========================================================
#
# Not really for export, but for aligning the binary data of BioSequences
# which don't have their binary data aligned. This is required for looping over
# two sequences integer by integer during mutation counting to do the
# bit parallel counting operations, but this may be useful in other instances to.
# Hence the iteration behaviour has been abstracted into this ShiftedInts iterator.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable ShiftedInts{A<:Alphabet}
    seq::BioSequence{A}
    finalIndex::BitIndex
    finalMask::UInt64
    nIntegers::Int64

    function ShiftedInts{A}(seq::BioSequence{A})
        bi = bitindex(seq, endof(seq))
        return new(seq,
                   bi,
                   mask(offset(bi) + bitsof(A)),
                   seq_data_len(A, length(seq)))
    end
end

@inline function Base.start(itr::ShiftedInts)
    return bitindex(itr.seq, 1), 1
end

@inline function Base.next(itr::ShiftedInts, state::Tuple{BitIndex, Int64})
    firstIdx = index(state[1])
    # Determine if the state has reached the final integer containing sequence data.
    firstIsFinal = firstIdx == index(itr.finalIndex)

    # Determine then if a second integer is needed from the data, and determine
    # if this second integer is the final integer containing sequence data.
    secondIdx = ifelse(firstIsFinal, index(state[1]), index(state[1] + 64))
    secondIsFinal = secondIdx == index(itr.finalIndex)

    # Get the first and second integers from the sequence data.
    firstInt = itr.seq.data[firstIdx]
    secondInt = itr.seq.data[secondIdx]

    # If either of the two integers is the final integer containing sequence data,
    # we need to mask them to make sure the data really does end.
    firstInt &= ifelse(firstIsFinal, itr.finalMask, 0xffffffffffffffff)
    secondInt &= ifelse(secondIsFinal, itr.finalMask, 0xffffffffffffffff)

    firstInt = firstInt >> offset(state[1])

    firstInt |= ifelse(firstIsFinal, firstInt, secondInt << (64 - offset(state[1])))

    return firstInt, (state[1] + 64, state[2] + 1)
end

@inline function Base.done(itr::ShiftedInts, state::Tuple{BitIndex, Int64})
    return state[2] > itr.nIntegers
end

@inline function Base.eltype(::Type{ShiftedInts})
    return UInt64
end

@inline function Base.length(itr::ShiftedInts)
    return itr.nIntegers
end
