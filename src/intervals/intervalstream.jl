
using Base.Collections: heappush!, heappop!

# IntervalStreams
# ---------------
#
# Often we'd rather avoid reading a large interval dataset into an
# IntervalCollection. We can still do efficient intersections if we can assume
# the data is sorted


type IntervalStreamIntersectIterator{S, T, SS, TS}
    a::SS
    b::TS
    seqname_isless::Function

    # buffered intervals from a
    a_buffer::StreamBuffer{Interval{S}}
    b_interval::Interval{T}

    function IntervalStreamIntersectIterator(a::SS,
                                             b::TS,
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
function intersect{SS <: IntervalStreamOrArray, TS <: IntervalStreamOrArray}(
        a::SS, b::TS, seqname_isless::Function=isless)
    S = metadatatype(a)
    T = metadatatype(b)
    return IntervalStreamIntersectIterator{S, T, SS, TS}(a, b, seqname_isless)
end


function first_intersection{S, T, SS, TS, U, V}(
                                        a::SS,
                                        b::TS,
                                        a_state::U,
                                        b_state::V,
                                        a_interval::Interval{S},
                                        b_interval::Interval{T},
                                        seqname_isless::Function)

    while !done(a, a_state) || !done(b, b_state)
        if precedes(a_interval, b_interval, seqname_isless) && !done(a, a_state)
            a_interval_next, a_state = next(a, a_state)
            if isless(a_interval_next, a_interval, seqname_isless)
                error("Interval streams must be sorted to perform intersection.")
            end
            a_interval = a_interval_next
        elseif precedes(b_interval, a_interval, seqname_isless) && !done(b, b_state)
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


function start{S, T, SS, TS}(it::IntervalStreamIntersectIterator{S, T, SS, TS})
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
function next{S, T, SS, TS, U, V}(it::IntervalStreamIntersectIterator{S, T, SS, TS},
                                  state::IntervalStreamIntersectIteratorState{U, V})
    value = (it.a_buffer[state.a_buffer_pos], it.b_interval)
    a_state = state.a_state
    b_state = state.b_state

    # another intersection in a_buffer?
    a_buffer_pos = state.a_buffer_pos + 1
    while a_buffer_pos <= length(it.a_buffer)
        if precedes(it.b_interval, it.a_buffer[a_buffer_pos], it.seqname_isless)
            a_buffer_pos = length(it.a_buffer) + 1
        elseif isoverlapping(it.b_interval, it.a_buffer[a_buffer_pos])
            break
        else
            a_buffer_pos += 1
        end
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
               precedes(b_interval, it.a_buffer[end], it.seqname_isless) ||
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
            #a_interval, b_interval, a_state, b_state = first_intersection(
                #it.a, it.b, a_state, b_state, a_interval, b_interval, it.seqname_isless)
            # See JuliaLang/julia#11181
            first_intersection_value = first_intersection(
                it.a, it.b, a_state, b_state, a_interval, b_interval, it.seqname_isless)
            a_interval = first_intersection_value[1]
            b_interval = first_intersection_value[2]
            a_state = first_intersection_value[3]
            b_state = first_intersection_value[4]

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


function done{S, T, SS, TS, U, V}(it::IntervalStreamIntersectIterator{S, T, SS, TS},
                                  state::IntervalStreamIntersectIteratorState{U, V})
    return state.a_buffer_pos > length(it.a_buffer)
end


# TODO: IntervalCollection(stream::IntervalStream) constructor.


# Helper function for coverage. Process remaining interval end points after
# all intervals have been read.
function coverage_process_lasts_heap!(cov::IntervalCollection{Uint32},
                                      current_coverage, coverage_seqname,
                                      coverage_first, lasts)
    while !isempty(lasts)
        pos = heappop!(lasts)
        if pos == coverage_first - 1
            current_coverage -= 1
        else
            @assert pos >= coverage_first
            push!(cov, Interval{Uint32}(coverage_seqname, coverage_first,
                                        pos, STRAND_BOTH, current_coverage))
            current_coverage -= 1
            coverage_first = pos + 1
        end
    end
    @assert current_coverage == 0
end

@doc """
Compute the coverage of a collection of intervals.

# Args
  * `intervals`: any IntervalStream

# Returns
An IntervalCollection that contains run-length encoded coverage data.

E.g. for intervals like
```
    [------]     [------------]
       [---------------]
```

this function would return a new set of disjoint intervals with annotated
coverage like:
```
    [1][-2-][-1-][--2--][--1--]
```
""" ->
function coverage(stream::IntervalStreamOrArray,
                  seqname_isless::Function=isless)
    cov = IntervalCollection{Uint32}()
    lasts = Int64[]

    stream_state = start(stream)
    if done(stream, stream_state)
        return cov
    end

    # TODO: How are we going to do this in a strand specific manner??
    # Just filter the stream?

    current_coverage = 0
    coverage_seqname = ""
    coverage_first = 0
    last_interval_first = 0
    interval, stream_state = next(stream, stream_state)
    while true
        if interval.seqname != coverage_seqname
            coverage_process_lasts_heap!(cov, current_coverage, coverage_seqname,
                                         coverage_first, lasts)
            if !(isempty(coverage_seqname) || seqname_isless(coverage_seqname, interval.seqname))
                error("Intervals must be sorted to compute coverage.")
            end

            coverage_seqname = interval.seqname
            current_coverage = 0
            coverage_first = 0
            last_interval_first = 0
        end

        if interval.first < last_interval_first
            error("Intervals must be sorted to compute coverage.")
        end

        if !isempty(lasts) && lasts[1] < first(interval)
            pos = heappop!(lasts)
            if first(interval) == pos + 1
                heappush!(lasts, last(interval))
                if done(stream, stream_state)
                    break
                end
                last_interval_first = first(interval)
                interval, stream_state = next(stream, stream_state)
            elseif pos == coverage_first - 1
                current_coverage -= 1
            else
                @assert pos >= coverage_first
                push!(cov, Interval{Uint32}(coverage_seqname, coverage_first,
                                            pos, STRAND_BOTH, current_coverage))
                current_coverage -= 1
                coverage_first = pos + 1
            end
        else
            if coverage_first == 0
                coverage_first = first(interval)
                current_coverage = 1
            elseif coverage_first == first(interval)
                current_coverage += 1
            else
                if current_coverage > 0
                    push!(cov, Interval{Uint32}(coverage_seqname, coverage_first,
                                                first(interval) - 1, STRAND_BOTH,
                                                current_coverage))
                end
                current_coverage += 1
                coverage_first = first(interval)
            end

            heappush!(lasts, last(interval))
            if done(stream, stream_state)
                break
            end
            last_interval_first = first(interval)
            interval, stream_state = next(stream, stream_state)
        end
    end

    coverage_process_lasts_heap!(cov, current_coverage, coverage_seqname,
                                 coverage_first, lasts)

    return cov
end


function coverage(ic::IntervalCollection)
    return coverage(ic, alphanum_isless)
end


