# Queue
# =====

type Queue{T}
    data::Vector{T}
    offset::Int
    first::Int
    last::Int
end

function (::Type{Queue{T}}){T}(bufsize::Integer=2^4)
    if bufsize ≤ 0
        throw(ArgumentError("buffer size must be positive"))
    elseif !ispow2(bufsize)
        throw(ArgumentError("buffer size must be a power of two"))
    end
    return Queue(Vector{T}(bufsize), 0, 1, 0)
end

function Base.eltype{T}(::Type{Queue{T}})
    return T
end

function Base.length(queue::Queue)
    return queue.last - queue.first + 1
end

function Base.push!{T}(queue::Queue{T}, elm::T)
    if length(queue.data) < length(queue) + 1
        index_first = dataindex(queue, queue.first)
        index_last = dataindex(queue, queue.last)
        index_end = endof(queue.data)
        # NOTE: resize factor must be a power of two
        resize!(queue.data, 2 * length(queue.data))
        @assert ispow2(length(queue.data))
        if !isempty(queue) && index_last < index_first
            # make the circular data linear
            copy!(queue.data, index_end + 1, queue.data, 1, index_last)
        end
        copy!(queue.data, 1, queue.data, index_first, length(queue))
        queue.offset = 0
    end
    queue.data[dataindex(queue, queue.last + 1)] = elm
    queue.last += 1
    return queue
end

function Base.shift!(queue::Queue)
    if isempty(queue)
        throw(ArgumentError("empty"))
    end
    elm = queue.data[dataindex(queue, queue.first)]
    queue.first += 1
    queue.offset += 1
    return elm
end

function Base.start(queue::Queue)
    return queue.first
end

function Base.done(queue::Queue, i)
    return i > queue.last
end

function Base.next(queue::Queue, i)
    return queue.data[dataindex(queue, i)], i + 1
end

function dataindex(queue::Queue, i::Integer)
    return ((i - queue.first + queue.offset) & (length(queue.data) - 1)) + 1
end
