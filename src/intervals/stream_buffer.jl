# Steam Buffer

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


function length{T}(buf::StreamBuffer{T})
    return buf.len
end


function isempty{T}(buf::StreamBuffer{T})
    return buf.last_idx == 0
end


function endof{T}(buf::StreamBuffer{T})
    return length(buf)
end


function push!{T}(buf::StreamBuffer{T}, x::T)
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


function shift!{T}(buf::StreamBuffer{T})
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


function getindex{T}(buf::StreamBuffer{T}, idx::Int)
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
