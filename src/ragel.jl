
module Ragel

using Compat
using Switch
import Base: push!, pop!, endof, append!, empty!, isempty, length, getindex,
             setindex!, (==), takebuf_string, read!, seek
using Bio: BufferedReader, fillbuffer!


# A simple buffer type, similar to IOBuffer but faster and with fewer features.
# It's essentially just a vector that can grow but is never allowed to
# reallocate to smaller size.

const INITIAL_BUF_SIZE = 1000000

type Buffer{T}
    data::Vector{T}
    pos::Int

    function Buffer(initial_size=10000)
        new(Array(T, initial_size), 1)
    end
end


function isempty(buf::Buffer)
    return buf.pos == 1
end


function length(buf::Buffer)
    return buf.pos - 1
end


function empty!(buf::Buffer)
    buf.pos = 1
end


function getindex(buf::Buffer, i::Integer)
    if !(1 <= i < buf.pos)
        throw(BoundsError())
    end
    @inbounds x = buf.data[i]
    return x
end


function setindex!{T}(buf::Buffer{T}, value::T, i::Integer)
    if !(1 <= i < buf.pos)
        throw(BoundsError())
    end
    @inbounds buf.data[i] = value
end


function ensureroom!(buf::Buffer, n::Integer)
    if buf.pos + n > length(buf.data)
        newsize = max(2 * length(buf.data), buf.pos + n)
        resize!(buf.data, newsize)
    end
end


function push!{T}(buf::Buffer{T}, c::T)
    ensureroom!(buf, 1)
    @inbounds buf.data[buf.pos] = c
    buf.pos += 1
end


function pop!(buf::Buffer)
    if buf.pos == 1
        throw(BoundsError())
    end
    @inbounds c = buf.data[buf.pos - 1]
    buf.pos -= 1
    return c
end


function endof(buf::Buffer)
    return buf.pos - 1
end


function append!{T}(buf::Buffer{T}, source::Vector{T}, start::Integer, stop::Integer)
    n = stop - start + 1
    ensureroom!(buf, n)
    copy!(buf.data, buf.pos, source, start, n)
    buf.pos += n
end


function (==){T}(a::Buffer{T}, b::Buffer{T})
    if a.pos != b.pos
        return false
    end
    for i in 1:a.pos-1
        if a.data[i] != b.data[i]
            return false
        end
    end
    return true
end


function (==){T}(a::Buffer{T}, b::String)
    if length(a) != length(b)
        return false
    end
    for (i, u) in enumerate(b)
        v = a.data[i]
        if u != v
            return false
        end
    end

    return true
end


function takebuf_string{T}(buf::Buffer{T})
    s = bytestring(buf.data[1:buf.pos-1])
    empty!(buf)
    return s
end


# Nuber of bytes to read at a time
const RAGEL_PARSER_INITIAL_BUF_SIZE = 1000000

# A type keeping track of a ragel-based parser's state.
type State
    reader::BufferedReader

    # true when all input has been parsed
    finished::Bool

    # Internal ragel state:
    p::Int  # index into the input stream (0-based)
    cs::Int # current DFA stae
    dfa_start_state::Int # Initial DFA state

    # Parser is responsible for updating this
    linenum::Int

    function State(cs, data::Vector{Uint8})
        return new(nothing, false, data, Buffer{Int}(16), false, 0,
                   length(data), cs, cs, 1)
    end

    function State(cs, input::IO, memory_map=false)
        if memory_map
            error("Parser must be given a file name in order to memory map.")
        end
        return new(BufferedReader(input, memory_map), false, 0, cs, cs, 1)
    end

    function State(cs, filename::String, memory_map=false)
        return new(BufferedReader(filename, memory_map), false, 0,  cs, cs, 1)
    end
end


# Get a State object from a parser. Parser implementations may want
# to define a more specific method.
function ragelstate(x)
    return x.state
end


# Macros for push and popping marks from within a ragel parser
macro mark!()
    quote
        $(esc(:state)).reader.mark = 1 + $(esc(:p))
    end
end

macro unmark!()
    quote
        m = $(esc(:state)).reader.mark
        $(esc(:state)).reader.mark = 0
        m
    end
end

macro position()
    quote
        1 + $(esc(:p))
    end
end

macro spanfrom(firstpos)
    quote
        $(esc(:state)).reader.buffer[$(esc(firstpos)):$(esc(:p))]
    end
end

macro asciistring_from_mark!()
    quote
        firstpos = Ragel.@unmark!
        ASCIIString($(esc(:state)).reader.buffer[firstpos:$(esc(:p))])
    end
end


# Parse an integer from an interval in the input buffer. This differs from
# `parse(Int, ...)` in that we don't have to copy or allocate anything, and it
# doesn't check that the characters are digits (we don't need to since this is
# already checked during parsing).
function parse_int64(buffer, firstpos, lastpos)
    x = @compat Int64(0)
    for i in firstpos:lastpos
        x = x * 10 + buffer[i] - (@compat UInt8('0'))
    end
    return x
end


macro int64_from_mark!()
    quote
        firstpos = @unmark!
        parse_int64($(esc(:state)).reader.buffer, firstpos, $(esc(:p)))
    end
end


# return the current character
macro char()
    quote
        convert(Char, $(esc(:state)).reader.buffer[$(esc(:p))+1])
    end
end


# Define a read! function wrapping ragel-generated parser.
#
# This macro handles some the dirty work of maintaining state, refilling
# buffers, etc.
#
# Args:
#
macro generate_read_fuction(machine_name, input_type, output_type, ragel_body, accept_body)
    start_state = esc(symbol(string(machine_name, "_start")))
    accept_state = esc(symbol(string(machine_name, "_first_final")))
    error_state = esc(symbol(string(machine_name, "_error")))

    # ragel needs these specific variable names so we have to escape them
    # throughout
    p = esc(:p)
    pe = esc(:pe)
    cs = esc(:cs)
    data = esc(:data)
    state = esc(:state)
    eof = esc(:eof)
    yield = esc(:yield)

    quote
        function $(esc(:advance!))(input::$(esc(input_type)))
            local $(esc(:input)) = input

            $(state) = ragelstate(input)
            if $(state).finished
                return false
            end

            $(p) = $(state).p
            $(pe) = $(state).reader.buffer_end
            $(cs) = $(state).cs
            $(data) = $(state).reader.buffer
            $(yield) = false

            # run the parser until all input is consumed or a match is found
            local $(eof) = $(pe) + 1
            while true
                if $(p) == $(pe)
                    $(state).p = $(p)
                    $(state).reader.buffer_end = $(pe)
                    nb = fillbuffer!($(state).reader)
                    $(p) = $(state).p
                    $(pe) = $(state).reader.buffer_end
                    if nb == 0
                        $(eof) = $(pe) # trigger ragel's eof handling
                    else
                        $(eof) = $(pe) + 1
                    end
                end

                $(esc(ragel_body))

                if $(cs) == $(error_state) && !$(yield)
                    error(string(
                        $("Error parsing $(machine_name) input on line "),
                        $(state).linenum))
                elseif $(yield)
                    if $(p) == $(pe)
                        $(state).p = $(p)
                        $(state).pe = $(pe)
                        fillbuffer!($(state).reader) == 0
                        $(p) = $(state).p
                        $(pe) = $(state).reader.buffer_end
                    end

                    break
                elseif $(p) == $(pe) == $(eof)
                    break
                end
            end

            if $(p) == $(pe) && $(cs) < $(accept_state)
                error($("Unexpected end of input while parsing $(machine_name)"))
            end

            $(state).p = $(p)
            $(state).cs = $(cs)

            if $(p) >= $(pe)
                $(state).finished = true
            end
            return true
        end
    end
end


end # module Ragel

