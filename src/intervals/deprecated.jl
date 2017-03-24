import Base.@deprecate, Base.@deprecate_binding
import Base.depwarn

# Bio.jl v0.4
# -----------

# TODO: any way to deprecate intersect?
#@deprecate Base.intersect(bb::BigBedReader, query::Interval)                       eachoverlap(bb, query)
#@deprecate Base.intersect{T}(a::IntervalCollection{T}, b::Interval)                eachoverlap(a, b)
#@deprecate Base.intersect{S,T}(a::IntervalCollection{S}, b::IntervalCollection{T}) eachoverlap(a, b)
#@deprecate Base.intersect(a::IntervalCollection, b::IntervalStreamOrArray)         eachoverlap(a, b)

"""
A type deriving `IntervalStream{T}` must be iterable and produce
Interval{T} objects in sorted order.
"""
abstract IntervalStream{T}
export IntervalStream

typealias IntervalStreamOrArray{T} Union{Vector{Interval{T}},IntervalStream{T},Bio.IO.AbstractReader}

# NOTE: Copied from src/intervals/streambuffer.jl
# Steam Buffer
# ============
#
# Stream buffer type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A buffer in which pushing to the back, popping from the front (shift),
and indexing is efficient.

It's constructed as a circular buffer that resizes when it runs out of space.
"""
type StreamBuffer{T}
    data::Vector{T}
    first_idx::Int
    last_idx::Int
    len::Int

    function StreamBuffer()
        return new(Array(T, 16), 1, 0, 0)
    end
end


function Base.length{T}(buf::StreamBuffer{T})
    return buf.len
end


function Base.isempty{T}(buf::StreamBuffer{T})
    return buf.last_idx == 0
end


function Base.endof{T}(buf::StreamBuffer{T})
    return length(buf)
end


function Base.push!{T}(buf::StreamBuffer{T}, x::T)
    if length(buf) == length(buf.data)
        newdata = Array(T, 2 * length(buf.data))
        len = length(buf)
        for i in 1:len
            newdata[i] = buf[i]
        end
        buf.data = newdata
        buf.first_idx = 1
        buf.last_idx = len
    end

    buf.last_idx += 1
    if buf.last_idx > length(buf.data)
        buf.last_idx -= length(buf.data)
    end
    @inbounds buf.data[buf.last_idx] = x
    buf.len += 1
    return x
end


function Base.shift!{T}(buf::StreamBuffer{T})
    if buf.last_idx == 0
        throw(BoundsError())
    end

    @inbounds x = buf.data[buf.first_idx]

    if buf.first_idx == buf.last_idx
        buf.last_idx = 0
        buf.first_idx = 1
    elseif buf.first_idx == length(buf.data)
        buf.first_idx = 1
    else
        buf.first_idx += 1
    end

    buf.len -= 1
    return x
end


function Base.getindex{T}(buf::StreamBuffer{T}, idx::Int)
    if idx > length(buf) || idx < 1
        throw(BoundsError())
    end

    idx = buf.first_idx + idx - 1
    if idx > length(buf.data)
        idx -= length(buf.data)
    end
    @inbounds x = buf.data[idx]
    return x
end


# NOTE: Copied from src/intervals/intervalstream.jl.
# IntervalStreams
# ===============
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

function Base.iteratorsize{S,T,SS,TS}(::Type{IntervalStreamIntersectIterator{S,T,SS,TS}})
    return Base.SizeUnknown()
end

immutable IntervalStreamIntersectIteratorState{U, V}
    a_buffer_pos::Int
    a_state::U
    b_state::V
end

"""
Intersect two `IntervalStreams` returning an iterator over all pairs of
intersecting intervals.
"""
function Base.intersect{SS<:IntervalStreamOrArray,TS<:IntervalStreamOrArray}(
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

function Base.start{S,T,SS,TS}(it::IntervalStreamIntersectIterator{S,T,SS,TS})
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
function Base.next{S,T,SS,TS,U,V}(it::IntervalStreamIntersectIterator{S,T,SS,TS},
                                  state::IntervalStreamIntersectIteratorState{U,V})
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
            if !isordered(it.a_buffer[end], a_interval, it.seqname_isless)
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
        if !isordered(b_interval, b_interval_next, it.seqname_isless)
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

function Base.done{S,T,SS,TS,U,V}(it::IntervalStreamIntersectIterator{S,T,SS,TS},
                                  state::IntervalStreamIntersectIteratorState{U,V})
    return state.a_buffer_pos > length(it.a_buffer)
end


# Bio.jl v0.3
# -----------

immutable BigWig <: Bio.IO.FileFormat end
immutable BigBed <: Bio.IO.FileFormat end
export BigWig, BigBed
