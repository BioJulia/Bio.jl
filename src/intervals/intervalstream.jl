

# IntervalStreams
# ---------------
#
# Often we'd rather avoid reading a large interval dataset into an
# IntervalCollection. We can still do efficient intersections if we can assume
# the data is sorted


@doc """
A type deriving `IntervalStream{T}` must be iterable and produce
Interval{T} objects in sorted order.
""" ->
abstract IntervalStream{T}


typealias IntervalStreamOrArray{T} Union(AbstractArray{Interval{T}}, IntervalStream{T})


type IntervalStreamIntersectIterator{S, T}
    a::IntervalStreamOrArray{S}
    b::IntervalStreamOrArray{T}
    seqname_isless::Function

    # buffered intervals from a
    a_buffer::StreamBuffer{Interval{S}}
    b_interval::Interval{T}

    function IntervalStreamIntersectIterator(a::IntervalStreamOrArray{S},
                                             b::IntervalStreamOrArray{T},
                                             seqname_isless::Function)
        return new(a, b, seqname_isless, StreamBuffer{Interval{S}}())
    end
end


immutable IntervalStreamIntersectIteratorState{U, V}
    a_buffer_pos::Int
    a_state::U
    b_state::V
end


@doc """
Intersect two `IntervalStreams` returning an iterator over all pairs of
intersecting intervals.
""" ->
function intersect{S, T}(a::IntervalStreamOrArray{S},
                         b::IntervalStreamOrArray{T},
                         seqname_isless::Function=isless)
    return IntervalStreamIntersectIterator{S, T}(a, b, seqname_isless)
end


function first_intersection{S, T, U, V}(a::IntervalStreamOrArray{S},
                                        b::IntervalStreamOrArray{T},
                                        a_state::U,
                                        b_state::V,
                                        a_interval::Interval{S},
                                        b_interval::Interval{T},
                                        seqname_isless::Function)

    while !done(a, a_state) || !done(b, b_state)
        if precedes(a_interval, b_interval) && !done(a, a_state)
            a_interval_next, a_state = next(a, a_state)
            if isless(a_interval_next, a_interval, seqname_isless)
                error("Interval streams must be sorted to perform intersection.")
            end
            a_interval = a_interval_next
        elseif precedes(b_interval, a_interval) && !done(b, b_state)
            b_interval_next, b_state = next(b, b_state)
            if isless(b_interval_next, b_interval, seqname_isless)
                error("Interval streams must be sorted to perform intersection.")
            end
            b_interval = b_interval_next
        else # intersection
            break
        end
    end

    return a_interval, b_interval, a_state, b_state
end


function start{S, T}(it::IntervalStreamIntersectIterator{S, T})
    a_state = start(it.a)
    b_state = start(it.b)

    if done(it.a, a_state) || done(it.b, b_state)
        return IntervalStreamIntersectIteratorState(a_state, b_state)
    end

    a_interval, a_state = next(it.a, a_state)
    b_interval, b_state = next(it.b, b_state)
    a_interval, b_interval, a_state, b_state = first_intersection(
        it.a, it.b, a_state, b_state, a_interval, b_interval, it.seqname_isless)

    if isoverlapping(a_interval, b_interval)
        push!(it.a_buffer, a_interval)
        it.b_interval = b_interval
    end

    return IntervalStreamIntersectIteratorState(1, a_state, b_state)
end


# Do we need a state if this is going to be stateful. State can maybe just
# be the current interval in b
function next{S, T, U, V}(it::IntervalStreamIntersectIterator{S, T},
                          state::IntervalStreamIntersectIteratorState{U, V})
    value = (it.a_buffer[state.a_buffer_pos], it.b_interval)
    a_state = state.a_state
    b_state = state.b_state

    # another intersection in a_buffer?
    a_buffer_pos = state.a_buffer_pos + 1
    while a_buffer_pos <= length(it.a_buffer) &&
          !isoverlapping(it.a_buffer[a_buffer_pos], it.b_interval)
        a_buffer_pos += 1
    end

    # look for an intersection for b_interval by reading more entries from a
    b_interval = it.b_interval
    if a_buffer_pos > length(it.a_buffer) && !done(it.a, a_state)
        while true
            a_interval, a_state = next(it.a, a_state)
            if isless(a_interval, it.a_buffer[end], it.seqname_isless)
                error("Interval streams must be sorted to perform intersection.")
            end
            push!(it.a_buffer, a_interval)

            if isoverlapping(b_interval, it.a_buffer[end]) ||
               precedes(b_interval, it.a_buffer[end]) ||
               done(it.a, a_state)
                break
            end
            a_buffer_pos += 1
        end

        if !isoverlapping(it.a_buffer[a_buffer_pos], b_interval)
            a_buffer_pos += 1
        end
    end

    # no more intersections against b_interval, try next interval from b
    while (a_buffer_pos > length(it.a_buffer) ||
           !isoverlapping(it.a_buffer[a_buffer_pos], b_interval)) && !done(it.b, b_state)

        b_interval_next, b_state = next(it.b, b_state)
        if isless(b_interval_next, b_interval, it.seqname_isless)
            error("Interval streams must be sorted to perform intersection.")
        end
        b_interval = b_interval_next

        # evict as much as possible from a_buffer
        a_interval = it.a_buffer[1]
        while !isempty(it.a_buffer) && precedes(it.a_buffer[1], b_interval, it.seqname_isless)
            a_interval = shift!(it.a_buffer)
        end

        if isempty(it.a_buffer)
            a_interval, b_interval, a_state, b_state = first_intersection(
                it.a, it.b, a_state, b_state, a_interval, b_interval, it.seqname_isless)
            it.b_interval = b_interval
            if isoverlapping(a_interval, b_interval)
                push!(it.a_buffer, a_interval)
                a_buffer_pos = 1
            end
            break
        elseif isoverlapping(it.a_buffer[1], b_interval)
            a_buffer_pos = 1
            it.b_interval = b_interval
        end
    end

    return value, IntervalStreamIntersectIteratorState{U, V}(
                        a_buffer_pos, a_state, b_state)
end


function done{S, T, U, V}(it::IntervalStreamIntersectIterator{S, T},
                          state::IntervalStreamIntersectIteratorState{U, V})
    return state.a_buffer_pos > length(it.a_buffer)
end


# TODO: IntervalCollection(stream::IntervalStream) constructor.
