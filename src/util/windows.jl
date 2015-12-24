
arrayOrStringOrSeq = Union{AbstractArray, AbstractString, Sequence}


"""
An iterator which moves across a container, as a sliding window.
"""
immutable EachWindowIterator{T <: arrayOrStringOrSeq}
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
function size(winitr::EachWindowIterator)
    return length(StepRange(winitr.width, winitr.step, length(winitr.data)))
end

"""
Calculate the number of elements in the container that are missed during iteration.
Typically because of the width and step size of the window: some elements near the end of the
container are missed if window sizes and step sizes are too large.
"""
function missed(winitr::EachWindowIterator)
    l = length(winitr.data)
    return l - StepRange(winitr.width, winitr.step, l)
end

"""
Convienient function for constructing an EachWindowIterator.

Accepts the same arguments as the EachWindowIterator constructor.
"""
function eachwindow{T <: arrayOrStringOrSeq}(data::T, width::Int, step::Int = 1)
    EachWindowIterator{T}(data, width, step)
end

@inline function start(it::EachWindowIterator)
    return 1
end

@inline function next(it::EachWindowIterator, state::Integer)
    i = state
    window = sub(it.data, i:i + it.width - 1)
    return window, i + it.step
end

@inline function done(it::EachWindowIterator, state::Integer)
    return (state + it.width - 1) > length(it.data)
end

function show(io::IO, it::EachWindowIterator)
    print(io, "Sliding window iterator.\nContainer size: ", length(it.data), " elements.",
    "\nWindow width: ", it.width, "\nStep distance: ", it.step,
    "\nTotal windows: ", size(it), "\nMissed elements: ",
    missed(it))
end
