# Iterate through every window of size K in vec
immutable EachWindowIterator{T}
    vec::Vector{T}
    K::Int
    step::Int
    npos::Int  # Number of windows, length(vec) - K + step

    EachWindowIterator(vec, K, step=1) = new(vec, K, step, length(vec) - K + 1)
end

immutable EachWindowIteratorState{T}
    i::Int  # start position of current window
end

function eachwindow{T}(K::Integer, vec::Vector{T}, step::Integer=1)
    @assert K >= 0 "K (window size) must be ≥ 0 in EachWindow"
    @assert step >= 1 "step must be ≥ 1"

    return EachWindowIterator{T}(vec, K, step)
end

@inline function start{T}(it::EachWindowIterator{T})
    i = 1
    npos = length(it.vec) - it.K + it.step
    return i
end

@inline function next{T}(it::EachWindowIterator{T}, state::Integer)
    i = state
    window = sub(it.vec, i:i+it.K - 1)
    return window, i + it.step
end

@inline function done{T}(it::EachWindowIterator{T}, state::Integer)
    return state > it.npos
end
