module Ragel

using BufferedStreams,
    Switch

using Bio:
    FileFormat,
    AbstractParser

import Base:
    push!,
    pop!,
    endof,
    append!,
    empty!,
    isempty,
    length,
    getindex,
    setindex!,
    (==),
    takebuf_string,
    read!,
    seek,
    open,
    start,
    next,
    done


# A type keeping track of a ragel-based parser's state.
type State{T <: BufferedInputStream}
    stream::T

    # Internal ragel state:
    p::Int  # index into the input stream (0-based)
    cs::Int # current DFA stae
    dfa_start_state::Int # Initial DFA state

    # Parser is responsible for updating this
    linenum::Int

    # true when all input has been consumed
    finished::Bool
end


"""
Construct a new ragel parser state.

### Args
  * `cs`: initial state
  * `input`: input stream
"""
function State{T <: BufferedInputStream}(cs::Int, input::T)
    if input.available == 0
        BufferedStreams.fillbuffer!(input)
    end

    return State{T}(input, 0, cs, cs, 1, false)
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


@inline function anchor!(state, p)
    state.stream.anchor = 1 + p
end


@inline function upanchor!(state)
    @assert state.stream.anchor != 0  "upanchor! called with no anchor set"
    anchor = state.stream.anchor
    state.stream.anchor = 0
    return anchor
end


# Parse an integer from an interval in the input buffer. This differs from
# `parse(Int, ...)` in that we don't have to copy or allocate anything, and it
# doesn't check that the characters are digits (we don't need to since this is
# already checked during parsing).
function parse_int64(buffer, firstpos, lastpos)
    x = Int64(0)
    for i in firstpos:lastpos
        x = x * 10 + buffer[i] - UInt8('0')
    end
    return x
end


# Macros that help make common parsing tasks moce succinct


macro anchor!()
    quote
        anchor!($(esc(:state)), $(esc(:p)))
    end
end


macro copy_from_anchor!(dest)
    quote
        firstpos = upanchor!($(esc(:state)))
        copy!($(esc(dest)), $(esc(:state)).stream.buffer, firstpos, $(esc(:p)))
    end
end


macro append_from_anchor!(dest)
    quote
        firstpos = upanchor!($(esc(:state)))
        append!($(esc(dest)), $(esc(:state)).stream.buffer, firstpos, $(esc(:p)))
    end
end


macro int64_from_anchor!()
    quote
        firstpos = upanchor!($(esc(:state)))
        parse_int64($(esc(:state)).stream.buffer, firstpos, $(esc(:p)))
    end
end


macro load_from_anchor!(T)
    quote
        firstpos = upanchor!($(esc(:state)))
        unsafe_load(convert(Ptr{$(T)},
            pointer($(esc(:state)).stream.buffer, firstpos)))
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
macro generate_read_fuction(machine_name, input_type, output_type, ragel_body)
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
    output = esc(:output)
    input = esc(:input)

    quote
        function $(esc(:(Base.read!)))(input::$(esc(input_type)),
                                     state::State, output::$(esc(output_type)))
            $(state) = state
            $(input) = input
            $(output) = output

            if $(state).finished
                return false
            end

            $(p) = $(state).p
            $(pe) = $(state).stream.available
            $(cs) = $(state).cs
            $(data) = $(state).stream.buffer
            $(yield) = false

            # run the parser until all input is consumed or a match is found
            $(eof) = $(pe) + 1
            @inbounds while true
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
            return $(yield)
        end

        function $(esc(:(Base.read!)))(input::$(esc(input_type)), output::$(esc(output_type)))
            # specialize on input.state
            $(esc(:read!))(input, input.state, output)
        end
    end
end


# Open functions for various sources
# ----------------------------------

function open{T <: FileFormat}(filename::AbstractString, ::Type{T}; args...)
    source = open(filename)
    stream = BufferedInputStream(source)
    open(stream, T; args...)
end


function open{T <: FileFormat}(filename::AbstractString, ::Type{T})
    source = open(filename)
    stream = BufferedInputStream(source)
    open(stream, T)
end


function open{T <: FileFormat}(source::Union{IO, Vector{UInt8}}, ::Type{T}; args...)
    open(BufferedInputStream(source), T; args...)
end


function open{T <: FileFormat}(source::Union(IO, Vector{UInt8}), ::Type{T})
    open(BufferedInputStream(source), T)
end


# Iterators for parsers
# ---------------------

function start{PT <: AbstractParser}(parser::PT)
    ET = eltype(PT)
    nextitem = ET()
    if read!(parser, nextitem)
        return Nullable{ET}(nextitem)
    else
        return Nullable{ET}()
    end
end


function next{ET}(parser::AbstractParser, nextitem_::Nullable{ET})
    nextitem = get(nextitem_)
    value = copy(nextitem)
    return (value,
        read!(parser, nextitem) ? Nullable{ET}(nextitem) : Nullable{ET}())
end


function done(parser::AbstractParser, nextitem::Nullable)
    return isnull(nextitem)
end


end # module Ragel
