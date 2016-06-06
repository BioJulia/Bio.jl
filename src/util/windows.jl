# Windows
# =======
#
# Sliding window iterator.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Windows

using Bio.Seq

using Compat

export eachwindow,
    EachWindowIterator,
    missed



typealias ArrayOrStringOrSeq Union{AbstractArray, AbstractString, BioSequence}

"""
An iterator which moves across a container, as a sliding window.

Constructor requires the data to move the window across, the width of the
window, and the step of the window.
"""
immutable EachWindowIterator{T <: ArrayOrStringOrSeq}
    "A reference to the collection to iterate over."
    data::T
    "The width of the sliding window."
    width::Int
    "How many positions the sliding window moves along."
    step::Int
    function EachWindowIterator{T}(data::T, width::Int, step::Int = 1)
        @assert width >= 1 "Window width must be ≥ 1."
        @assert step >= 1 "step must be ≥ 1."
        @assert width <= length(data) "The window size cannot be greater than number of data elements."
        return new(data, width, step)
    end
end


"""
Calculate the number of windows that will result from iterating across the container.

Accepts one variable of type EachWindowIterator.
"""
function Base.size(winitr::EachWindowIterator)
    return length(StepRange(winitr.width, winitr.step, length(winitr.data)))
end

Base.length(winitr::EachWindowIterator) = size(winitr)

"""
Calculate the number of elements in the container that are missed during iteration.
Typically because of the width and step size of the window: some elements near the end of the
container are missed if window sizes and step sizes are too large.
"""
function missed(winitr::EachWindowIterator)
    l = length(winitr.data)
    return l - StepRange(winitr.width, winitr.step, l).stop
end

"""
A sliding window iterator is considered equal to another
sliding window iterator if the data is considered equal, and the
width and step of the iterator is also equivalent.
"""
@compat function Base.:(==)(x::EachWindowIterator, y::EachWindowIterator)
    return (x.data == y.data) && (x.width == y.width) && (x.step == y.step)
end

"""
Convienient function for constructing an EachWindowIterator.

Accepts the same arguments as the EachWindowIterator constructor.
"""
function eachwindow{T <: ArrayOrStringOrSeq}(data::T, width::Int, step::Int = 1)
    EachWindowIterator{T}(data, width, step)
end

@inline function Base.start(it::EachWindowIterator)
    return 1
end

@inline function Base.next(it::EachWindowIterator, state::Integer)
    i = state
    window = sub(it.data, i:i + it.width - 1)
    return (i, window), i + it.step
end

# Extra next method to account for fact than strings don't have sub method
# like arrays and Nucleotide sequences do.
# The operation is an indexing operation, rather than a "proper"
# substring or subarray operation.
@inline function Base.next{T <: AbstractString}(it::EachWindowIterator{T}, state::Integer)
    i = state
    window = it.data[i:i + it.width - 1]
    return (i, window), i + it.step
end

@inline function Base.done(it::EachWindowIterator, state::Integer)
    return (state + it.width - 1) > length(it.data)
end

function Base.show(io::IO, it::EachWindowIterator)
    print(io, "Sliding window iterator.\nContainer size: ", length(it.data), " elements.",
    "\nWindow width: ", it.width, "\nStep distance: ", it.step,
    "\nTotal windows: ", size(it), "\nMissed elements: ",
    missed(it))
end

end
