# Ragel
# =====
#
# Utilities for the Ragel state machine compiler.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Ragel

export tryread!

using BufferedStreams
import Bio.IO: FileFormat, AbstractReader, stream

# A type keeping track of a ragel-based parser's state.
type State{T<:BufferedInputStream}
    stream::T      # input stream
    cs::Int        # current DFA state of Ragel
    linenum::Int   # line number: parser is responsible for updating this
    finished::Bool # true if finished (regardless of where in the stream we are)
end

function State(initstate::Int, input::BufferedInputStream)
    return State(input, initstate, 1, false)
end

@inline function anchor!(state, p)
    state.stream.anchor = 1 + p
end

@inline function upanchor!(state)
    @assert state.stream.anchor != 0 "upanchor! called with no anchor set"
    anchor = state.stream.anchor
    state.stream.anchor = 0
    return anchor
end

# Parse an integer from an interval in the input buffer. This differs from
# `parse(Int, ...)` in that we don't have to copy or allocate anything, and it
# doesn't check that the characters are digits (we don't need to since this is
# already checked during parsing).
function parse_int64(buffer, firstpos, lastpos)
    if buffer[firstpos] == UInt8('-')
        sign = -1
        firstpos += 1
    else
        sign = +1
    end
    x = Int64(0)
    for i in firstpos:lastpos
        x = x * 10 + buffer[i] - UInt8('0')
    end
    return sign * x
end


# Actions
# -------

# Macros that help make common parsing tasks more succinct

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

macro float64_from_anchor!()
    quote
        firstpos = upanchor!($(esc(:state))) - 1
        get(ccall(
            :jl_try_substrtod,
            Nullable{Float64},
            (Ptr{UInt8}, Csize_t, Csize_t),
            $(esc(:state)).stream.buffer, firstpos, $(esc(:p)) - firstpos
        ))
    end
end

macro ascii_from_anchor!()
    quote
        firstpos = upanchor!($(esc(:state)))
        n = $(esc(:p)) - firstpos + 1
        dst = Vector{UInt8}(n)
        copy!(dst, 1, $(esc(:state)).stream.buffer, firstpos, n)
        String(dst)
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

# advance `p`, set the current state to the target state `ts` and escape from
# scanning input
macro yield(ts)
    esc(quote
        p += 1
        cs = $ts
        @goto yield
    end)
end


# read/read! function generators
# ------------------------------

# Define a read!/read function wrapping ragel-generated parser.  These macros
# handle some the dirty work of maintaining state, refilling buffers, etc.
macro generate_read!_function(machine_name, input_type, output_type, ragel_body)
    accept_state = Symbol(machine_name, "_first_final")
    error_state = Symbol(machine_name, "_error")

    return esc(quote
        function Base.read!(input::$(input_type), output::$(output_type))
            state = input.state

            if state.finished
                throw(EOFError())
            end

            # restore state variables used by Ragel (`p` is 0-based)
            cs = state.cs
            p = state.stream.position - 1
            pe = state.stream.available
            eof = -1
            data = state.stream.buffer

            if !isopen(state.stream)
                error("input stream is already closed")
            end

            if p ≥ pe
                # fill data buffer
                if BufferedStreams.fillbuffer!(state.stream) == 0
                    state.finished = true
                    throw(EOFError())
                end
                p = state.stream.position - 1
                pe = state.stream.available
            end

            # run the parser until all input is consumed or a match is found
            while true
                $(ragel_body)

                state.cs = cs
                state.stream.position = p + 1

                if cs == $(error_state)
                    error("error parsing ", $(machine_name),
                          " input on line ", state.linenum)
                elseif p == eof || state.finished  # exhausted all input data
                    if cs >= $(accept_state)
                        state.finished = true
                        throw(EOFError())
                    else
                        error("unexpected end of input while parsing ", $(machine_name))
                    end
                elseif p == pe  # exhausted filled input data
                    # refill data buffer
                    hits_eof = BufferedStreams.fillbuffer!(state.stream) == 0
                    p = state.stream.position - 1
                    pe = state.stream.available
                    if hits_eof
                        eof = pe
                    end
                else
                    # unreachable here
                    @assert false
                end
            end

            # `@yield` jumps here
            @label yield

            # save the current state and the reading position
            state.cs = cs
            state.stream.position = p + 1

            return output
        end
    end)
end

function Base.eof(input::AbstractReader)
    return input.state.finished || eof(stream(input))
end

function Base.read(input::AbstractReader)
    return read!(input, eltype(input)())
end

"""
    tryread!(reader::AbstractReader, output)

Try to read the next element into `output` from `reader`.

The result is wrapped in `Nullable` and will be null if no entry is available.
"""
function tryread!(reader::AbstractReader, output)
    T = eltype(reader)
    try
        read!(reader, output)
        return Nullable{T}(output)
    catch ex
        if isa(ex, EOFError)
            return Nullable{T}()
        end
        rethrow()
    end
end


# Iterator
# --------

function Base.start(reader::AbstractReader)
    T = eltype(reader)
    nextone = T()
    if isnull(tryread!(reader, nextone))
        return Nullable{T}()
    else
        return Nullable{T}(nextone)
    end
end

Base.done(reader::AbstractReader, nextone) = isnull(nextone)

function Base.next(reader::AbstractReader, nextone)
    item = get(nextone)
    ret = copy(item)
    if isnull(tryread!(reader, item))
        return ret, Nullable{eltype(reader)}()
    else
        return ret, Nullable(item)
    end
end

end # module Ragel
