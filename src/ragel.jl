
module Ragel

using Compat
using BufferedStreams
using Switch
import Base: push!, pop!, endof, append!, empty!, isempty, length, getindex,
             setindex!, (==), takebuf_string, read!, seek


# A type keeping track of a ragel-based parser's state.
type State
    stream::BufferedInputStream

    # Internal ragel state:
    p::Int  # index into the input stream (0-based)
    cs::Int # current DFA stae
    dfa_start_state::Int # Initial DFA state

    # Parser is responsible for updating this
    linenum::Int

    # true when all input has been consumed
    finished::Bool

    function State(cs, data::Vector{Uint8}, memory_map::Bool=false, len::Integer=length(data))
        if memory_map
            error("Parser must be given a file name in order to memory map.")
        end
        return new(BufferedInputStream(data), 0, cs, cs, 1, false)
    end

    function State(cs, input::IO, memory_map=false)
        if memory_map
            error("Parser must be given a file name in order to memory map.")
        end
        stream = BufferedInputStream(input)
        BufferedStreams.fillbuffer!(stream)
        return new(stream, 0, cs, cs, 1, false)
    end

    function State(cs, filename::String, memory_map=false)
        if memory_map
            input = Mmap.mmap(open(filename), Vector{Uint8}, (filesize(filename),))
        else
            input = open(filename)
        end
        stream = BufferedInputStream(input)
        BufferedStreams.fillbuffer!(stream)
        return new(stream, 0, cs, cs, 1, false)
    end
end


function BufferedStreams.fillbuffer!(state::State)
    old_available = state.stream.available
    nb = BufferedStreams.fillbuffer!(state.stream)
    state.p = state.p + state.stream.available - old_available - nb
    return nb
end


# Get a State object from a parser. Parser implementations may want
# to define a more specific method.
function ragelstate(x)
    return x.state
end


# Macros for push and popping anchors from within a ragel parser
macro anchor!()
    quote
        $(esc(:state)).stream.anchor = 1 + $(esc(:p))
    end
end

macro upanchor!()
    quote
        @assert $(esc(:state)).stream.anchor != 0  "upanchor! called with no anchor set"
        a = $(esc(:state)).stream.anchor
        $(esc(:state)).stream.anchor = 0
        a
    end
end

macro position()
    quote
        1 + $(esc(:p))
    end
end

macro spanfrom(firstpos)
    quote
        $(esc(:state)).stream.buffer[$(esc(firstpos)):$(esc(:p))]
    end
end

macro asciistring_from_anchor!()
    quote
        firstpos = Ragel.@upanchor!
        ASCIIString($(esc(:state)).stream.buffer[firstpos:$(esc(:p))])
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


macro int64_from_anchor!()
    quote
        firstpos = @upanchor!
        parse_int64($(esc(:state)).stream.buffer, firstpos, $(esc(:p)))
    end
end


# return the current character
macro char()
    quote
        convert(Char, $(esc(:state)).stream.buffer[$(esc(:p))+1])
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
            $(pe) = $(state).stream.available
            $(cs) = $(state).cs
            $(data) = $(state).stream.buffer
            $(yield) = false

            # run the parser until all input is consumed or a match is found
            local $(eof) = $(pe) + 1
            while true
                if $(p) == $(pe)
                    $(state).p = $(p)
                    $(state).stream.available = $(pe)
                    nb = fillbuffer!($(state))
                    $(p) = $(state).p
                    $(pe) = $(state).stream.available
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
                        $(state).stream.available = $(pe)
                        fillbuffer!($(state))
                        $(p) = $(state).p
                        $(pe) = $(state).stream.available
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
            $(state).stream.available = $(pe)
            $(state).cs = $(cs)

            if $(p) >= $(pe)
                $(state).finished = true
            end
            return true
        end
    end
end


end # module Ragel

