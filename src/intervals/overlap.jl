# Overlap Iterator
# ================

immutable OverlapIterator{Sa,Sb,F}
    intervals_a::Sa
    intervals_b::Sb
    isless::F
end

function Base.eltype{Sa,Sb,F}(::Type{OverlapIterator{Sa,Sb,F}})
    return Tuple{Interval{metadatatype(Sa)},Interval{metadatatype(Sb)}}
end

function Base.iteratorsize{Sa,Sb,F}(::Type{OverlapIterator{Sa,Sb,F}})
    return Base.SizeUnknown()
end

"""
    eachoverlap(intervals_a, intervals_b, [seqname_isless=Base.isless])

Create an iterator of overlapping intervals between `intervals_a` and `intervals_b`.

This function assumes elements of `intervals_a` and `intervals_b` are sorted by
its sequence name and left position.  If the element type is not a subtype of
`Bio.Intervals.Interval`, elements are converted to `Interval` objects.

The third optional argument is a function that defines the order of sequence
names. The default function is `Base.isless`, which is the lexicographical
order.
"""
function eachoverlap(intervals_a, intervals_b, seqname_isless=Base.isless)
    return OverlapIterator(intervals_a, intervals_b, seqname_isless)
end

type OverlapIteratorState{Sa,Sb,Ta,Tb}
    state_a::Sa
    state_b::Sb
    queue::Queue{Interval{Tb}}
    state_queue::Int
    done::Bool
    interval_a::Interval{Ta}
    interval_b::Interval{Tb}

    function OverlapIteratorState(state_a, state_b)
        queue = Queue{Interval{Tb}}()
        return new(state_a, state_b, queue, start(queue), false)
    end
end

function Base.start(iter::OverlapIterator)
    state_a = start(iter.intervals_a)
    state_b = start(iter.intervals_b)
    Sa = typeof(state_a)
    Sb = typeof(state_b)
    Ta = metadatatype(iter.intervals_a)
    Tb = metadatatype(iter.intervals_b)
    state = OverlapIteratorState{Sa,Sb,Ta,Tb}(state_a, state_b)
    if done(iter.intervals_a, state_a) || done(iter.intervals_b, state_b)
        state.done = true
    else
        state.interval_a, state.state_a = next(iter.intervals_a, state_a)
        state.interval_b, state.state_b = next(iter.intervals_b, state_b)
        push!(state.queue, state.interval_b)
    end
    return advance!(state, iter)
end

function Base.done(iter::OverlapIterator, state)
    return state.done
end

function Base.next(iter::OverlapIterator, state)
    interval_a = state.interval_a
    interval_b, state.state_queue = next(state.queue, state.state_queue)
    return (interval_a, interval_b), advance!(state, iter)
end

function advance!(state::OverlapIteratorState, iter::OverlapIterator)
    # queue a new interval from intervals_b if possible
    function queue!()
        if done(iter.intervals_b, state.state_b)
            if !done(iter.intervals_a, state.state_a)
                elem_a, state.state_a = next(iter.intervals_a, state.state_a)
                next_interval_a = convert(typeof(state.interval_a), elem_a)
                check_ordered(state.interval_a, next_interval_a)
                state.interval_a = next_interval_a
                state.state_queue = start(state.queue)
            end
        else
            elem_b, state.state_b = next(iter.intervals_b, state.state_b)
            next_interval_b = convert(typeof(state.interval_b), elem_b)
            check_ordered(state.interval_b, next_interval_b)
            state.interval_b = next_interval_b
            push!(state.queue, next_interval_b)
        end
        return nothing
    end

    # check i1 and i2 are ordered
    function check_ordered(i1, i2)
        if !Bio.Intervals.isordered(i1, i2, iter.isless)
            error("intervals are not sorted")
        end
        return nothing
    end

    # find overlapping intervals if any and return there
    if state.done
        return state
    elseif done(state.queue, state.state_queue)
        queue!()
    end
    while !done(state.queue, state.state_queue)
        interval_b, next_state_queue = next(state.queue, state.state_queue)
        c = compare_overlap(state.interval_a, interval_b, iter.isless)
        if c < 0
            if done(iter.intervals_a, state.state_a)
                break
            else
                elem_a, state.state_a = next(iter.intervals_a, state.state_a)
                next_interval_a = convert(typeof(state.interval_a), elem_a)
                check_ordered(state.interval_a, next_interval_a)
                state.interval_a = next_interval_a
                state.state_queue = start(state.queue)
            end
        elseif c == 0
            return state
        else
            if interval_b === first(state.queue)
                shift!(state.queue)
            end
            state.state_queue = next_state_queue
            if done(state.queue, state.state_queue)
                queue!()
            end
        end
    end

    # no overlapping intervals found
    state.done = true
    return state
end

# Return:
#   -1 when `i1` precedes `i2`,
#   0 when `i1` overlaps with `i2`, and
#   +1 when `i1` follows `i2`.
function compare_overlap(i1::Interval, i2::Interval, isless::Function)
    if isless(i1.seqname, i2.seqname)::Bool
        return -1
    elseif isless(i2.seqname, i1.seqname)::Bool
        return +1
    else  # i1.seqname == i2.seqname
        if i1.last < i2.first
            return -1
        elseif i1.first > i2.last
            return +1
        else
            return 0
        end
    end
end

# Faster comparison for `Base.isless`.  Note that `Base.isless` must be
# consistent wtih `Base.cmp` to work correctly.
function compare_overlap(i1::Interval, i2::Interval, ::typeof(Base.isless))
    c = cmp(i1.seqname, i2.seqname)
    if c != 0
        return c
    end
    if i1.last < i2.first
        return -1
    elseif i1.first > i2.last
        return +1
    else
        return 0
    end
end
