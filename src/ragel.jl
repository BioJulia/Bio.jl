
module Ragel

using Compat
using Switch
import Base.FS
import Base: push!, pop!, endof, append!, empty!, isempty, length, getindex,
             setindex!, (==), takebuf_string, read!


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


function ensureroom!(buf::Buffer, n::Int)
    if buf.pos + n > length(buf.data)
        resize!(buf.data, 2 * length(buf.data))
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


function append!{T}(buf::Buffer{T}, source::Vector{T}, start::Int, stop::Int)
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
    # Input source. If nothing, it's assumed that buffer contains
    # all the input.
    input::Union(IO, FS.File, Nothing)

    # Should we close input when finished
    input_owned::Bool

    # Buffered data.
    buffer::Vector{Uint8}

    # Indexes into buffer that will be updated and stay valid upon refill.
    marks::Buffer{Int}

    # True when parsing has completed
    finished::Bool

    # Internal ragel state:
    p::Int  # index into the input stream (0-based)
    pe::Int # index after the last symbol in the input stream (0-based)
    cs::Int # current DFA stae

    # Parser is responsible for updating this
    linenum::Int

    function State(cs, data::Vector{Uint8})
        return new(nothing, false, data, Buffer{Int}(16), false, 0, length(data), cs, 1)
    end

    function State(cs, input::IO)
        return new(input, false, Array(Uint8, RAGEL_PARSER_INITIAL_BUF_SIZE),
                   Buffer{Int}(16), false, 0, 0, cs, 1)
    end

    function State(cs, filename::String, memory_map=false)
        if memory_map
            data = mmap_array(Uint8, (filesize(filename),), open(filename))
            return new(nothing, false, data, Buffer{Int}(16), false, 0, length(data), cs, 1)
        else
            file = FS.open(filename, FS.JL_O_RDONLY)
            return new(file, true, Array(Uint8, RAGEL_PARSER_INITIAL_BUF_SIZE),
                       Buffer{Int}(16), false, 0, 0, cs, 1)
        end
    end
end


# Get a State object from a parser. Parser implementations may want
# to define a more specific method.
function ragelstate(x)
    return x.state
end


# Macros for push and popping marks from within a ragel parser
macro pushmark!()
    quote
        push!($(esc(:state)).marks, 1 + $(esc(:p)))
    end
end

macro popmark!()
    quote
        pop!($(esc(:state)).marks)
    end
end

macro peekmark!()
    quote
        $(esc(:state)).marks[end]
    end
end

macro position()
    quote
        1 + $(esc(:p))
    end
end

macro spanfrom(firstpos)
    quote
        $(esc(:state)).buffer[$(esc(firstpos)):$(esc(:p))]
    end
end

macro asciistring_from_mark!()
    quote
        firstpos = Ragel.@popmark!
        ASCIIString($(esc(:state)).buffer[firstpos:$(esc(:p))])
    end
end

macro bytestring_from_mark!()
    quote
        firstpos = @popmark!
        len = $(esc(:p)) - firstpos + 1
        bytestring(pointer($(esc(:state)).buffer, firstpos), len)
    end
end


# Parse an integer from an interval in the input buffer. This differs from
# `parse(Int, ...)` in that we don't have to copy or allocate anything, and it
# doesn't check that the characters are digits (we don't need to since this is
# already checked during parsing).
function parse_int64(buffer, firstpos, lastpos)
    x = @compat Int64(0)
    for i in lastpos:-1:firstpos
        x = x * 10 + buffer[i] - '0'
    end
    return x
end


macro int64_from_mark!()
    quote
        firstpos = @popmark!
        parse_int64($(esc(:state)).buffer, firstpos, $(esc(:p)))
    end
end


# return the current character
macro char()
    quote
        convert(Char, $(esc(:state)).buffer[$(esc(:p))+1])
    end
end


# Refill a buffer, keeping some portion of it.
#
# The range buffer[parser.marks[1]:parser.pe] will be shifted to the beginning of
# the buffer, then the rest of the buffer will be filled by reading bytes from
# the input.
#
# This is useful if a parser is the middle of matching something when the buffer
# runs out
#
# Args:
#   parser: Current ragel parser state.
#
# Modifies:
#   parser.buffer is refilled, and indexes are updated to the correct positions
#   in the refilled buffer
#
# Returns:
#   Number of bytes read
#
function fillbuffer!(parser::State)
    if parser.input === nothing
        return 0
    end

    buflen = length(parser.buffer)
    keeplen = 0
    first_mark = 0
    if !isempty(parser.marks)
        first_mark = parser.marks[1]
        keeplen = parser.pe - first_mark + 1
        if keeplen == buflen
            buflen = 2 * buflen
            resize!(parser.buffer, buflen)
        end
        copy!(parser.buffer, 1, parser.buffer, first_mark, keeplen)
    end

    nb = readchunk!(parser.input, parser.buffer, keeplen + 1, buflen)

    parser.p = keeplen
    parser.pe = keeplen + nb
    for i in 1:length(parser.marks)
        parser.marks[i] -= first_mark - 1
    end

    if nb == 0 && parser.input_owned
        close(parser.input)
    end

    return nb
end


function readchunk!(source::IO, dest::Vector{Uint8}, dest_start::Int, dest_stop::Int)
    i = dest_start
    while i <= dest_stop && !eof(source)
        @inbounds dest[i] = read(source, Uint8)
        i += 1
    end
    return i - dest_start
end


function readchunk!(source::FS.File, dest::Vector{Uint8}, dest_start::Int, dest_stop::Int)
    len = dest_stop - dest_start  + 1
    return ccall(:jl_fs_read, Int32, (Int32, Ptr{Void}, Csize_t),
                 source.handle, pointer(dest, dest_start), len)
end


function readchunk!(source::IOStream, dest::Vector{Uint8}, dest_start::Int, dest_stop::Int)
    len = dest_stop - dest_start  + 1
    return ccall(:ios_readall, Uint, (Ptr{Void}, Ptr{Void}, Uint), source.ios,
                 pointer(dest, dest_start), len)
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

    quote
        function $(esc(:advance!))(input::$(esc(input_type)))
            # TODO: is there a more idiomatic way to do this?
            local $(esc(:input)) = input

            $(state) = ragelstate(input)
            if $(state).finished
                return false
            end

            $(p) = $(state).p
            $(pe) = $(state).pe
            $(cs) = $(state).cs
            $(data) = $(state).buffer
            $(esc(:yield)) = false

            # run the parser until all input is consumed or a match is found
            while true
                if $(p) == $(pe)
                    if fillbuffer!($(state)) == 0
                        break
                    end

                    $(p) = $(state).p
                    $(pe) = $(state).pe
                end

                $(esc(ragel_body))

                if $(cs) == $(error_state)
                    error(string(
                        $("Error parsing $(machine_name) input on line "),
                        $(state).linenum))
                elseif $(esc(:yield))
                    if $(p) == $(pe)
                        fillbuffer!($(state)) == 0
                        $(p) = $(state).p
                        $(pe) = $(state).pe
                    end

                    break
                end
            end

            if $(p) == $(pe) && $(cs) < $(accept_state)
                error($("Unexpected end of input while parsing $(machine_name)"))
            end

            $(state).p = $(p)
            $(state).pe = $(pe)
            $(state).cs = $(cs)

            if $(p) == $(pe)
                $(state).finished = true
            end
            return true
        end

        function read!(input::$(esc(input_type)), output::$(esc(output_type)))
            local $(esc(:input)) = input
            local $(esc(:output)) = output
            if advance!(input)
                $(esc(accept_body))
                return true
            else
                return false
            end
        end
    end
end


end # module Ragel

