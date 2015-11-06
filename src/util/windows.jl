
"""
An iterator which moves across a container, as a sliding window.
"""
immutable EachWindowIterator{T <: AbstractArray}
    "A reference to the collection to iterate over."
    data::T
    "The width of the sliding window."
    width::Int
    "How many positions the sliding window moves along."
    step::Int
    "The number of windows that will result from iterating across the container."
    nwin::Int
    "The number of elements in the container that are missed. Typically because
    of the width and step size of the window: some elements near the end of the
    container are missed if window sizes and step sizes are too large."
    nmissed::Int
end

function EachWindowIterator{T <: AbstractArray}(data::T, width::Int, step::Int = 1)
    @assert width >= 1 "Window width must be ≥ 1."
    @assert step >= 1 "step must be ≥ 1."
    @assert width < length(data) "The window size cannot be greater than number of data elements."
    r = StepRange(width, step, length(data))
    return EachWindowIterator(data, width, step, length(r), length(data) - last(r))
end

"""
Convienient function for constructing an EachWindowIterator.

Accepts the same arguments as the EachWindowIterator constructor.
"""
eachwindow{T <: AbstractArray}(data::T, width::Integer, step::Integer = 1) = EachWindowIterator(data, width, step)


@inline function start(it::EachWindowIterator)
    return 1
end

@inline function next{T <: AbstractArray}(it::EachWindowIterator{T}, state::Integer)
    i = state
    window = sub(it.data, i:i + it.width - 1)
    return window, i + it.step
end

@inline function done{T <: AbstractArray}(it::EachWindowIterator{T}, state::Integer)
    return (state + it.width - 1) > length(it.data)
end

function show{T <: AbstractArray}(io::IO, it::EachWindowIterator{T})
    print(io, "Sliding window iterator.\nContainer size: ", length(it.data), " elements.",
    "\nWindow width: ", it.width, "\nStep distance: ", it.step,
    "\nTotal windows: ", it.nwin, "\nMissed elements: ",
    it.nmissed)
end
