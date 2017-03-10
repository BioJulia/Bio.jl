# Reader Helper
# =============
#
# Utilities to generate file readers in Bio.jl.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module ReaderHelper

import Automa
import BufferedStreams

@inline function anchor!(stream::BufferedStreams.BufferedInputStream, p)
    stream.anchor = p
    stream.immobilized = true
    return stream
end

@inline function upanchor!(stream::BufferedStreams.BufferedInputStream)
    @assert stream.anchor != 0 "upanchor! called with no anchor set"
    anchor = stream.anchor
    stream.anchor = 0
    stream.immobilized = false
    return anchor
end

function ensure_margin!(stream::BufferedStreams.BufferedInputStream)
    if stream.position * 20 > length(stream.buffer) * 19
        BufferedStreams.shiftdata!(stream)
    end
    return nothing
end

function resize_and_copy!(dst::Vector{UInt8}, src::Vector{UInt8}, r::UnitRange{Int})
    rlen = length(r)
    if length(dst) != rlen
        resize!(dst, rlen)
    end
    copy!(dst, 1, src, first(r), rlen)
    return dst
end

function generate_index_function(record_type, machine, actions)
    quote
        function index!(record::$(record_type))
            data = record.data
            p = 1
            p_end = p_eof = sizeof(data)
            offset = mark = mark1 = mark2 = 0
            initialize!(record)
            cs = $(machine.start_state)
            $(Automa.generate_exec_code(machine, actions=actions, code=:goto, check=false))
            if cs != 0
                throw(ArgumentError(string("failed to index ", $(record_type), " ~>", repr(String(data[p:min(p+7,p_end)])))))
            end
            @assert isfilled(record)
            return record
        end
    end
end

function generate_readheader_function(reader_type, metainfo_type, machine, actions, finish_code=:())
    quote
        function readheader!(reader::$(reader_type))
            _readheader!(reader, reader.state)
        end

        function _readheader!(reader::$(reader_type), state::Bio.Ragel.State)
            stream = state.stream
            Bio.ReaderHelper.ensure_margin!(stream)
            cs = state.cs
            linenum = state.linenum
            data = stream.buffer
            p = stream.position
            p_end = stream.available
            p_eof = -1
            offset = mark1 = mark2 = 0
            finish_header = false
            record = $(metainfo_type)()

            while true
                $(Automa.generate_exec_code(machine, actions=actions, code=:table))

                state.cs = cs
                state.finished = cs == 0
                state.linenum = linenum
                stream.position = p

                if cs < 0
                    error("$($(reader_type)) file format error on line ", linenum)
                elseif finish_header
                    $(finish_code)
                    break
                elseif p > p_eof ≥ 0
                    error("incomplete $($(reader_type)) input on line ", linenum)
                else
                    hits_eof = BufferedStreams.fillbuffer!(stream) == 0
                    p = stream.position
                    p_end = stream.available
                    if hits_eof
                        p_eof = p_end
                    end
                end
            end
        end
    end
end

function generate_read_function(reader_type, machine, actions)
    quote
        function Base.read!(reader::$(reader_type), record::eltype($(reader_type)))::eltype($(reader_type))
            return _read!(reader, reader.state, record)
        end

        function _read!(reader::$(reader_type), state::Bio.Ragel.State, record::eltype($(reader_type)))
            stream = state.stream
            Bio.ReaderHelper.ensure_margin!(stream)
            initialize!(record)
            cs = state.cs
            linenum = state.linenum
            data = stream.buffer
            p = stream.position
            p_end = stream.available
            p_eof = -1
            offset = mark = 0
            found_record = false

            if state.finished
                throw(EOFError())
            end

            while true
                $(Automa.generate_exec_code(machine, actions=actions, code=:goto, check=false))

                state.cs = cs
                state.finished |= cs == 0
                state.linenum = linenum
                stream.position = p

                if cs < 0
                    error($(reader_type), " file format error on line ", linenum, " ~>", repr(String(data[p:min(p+7,p_end)])))
                elseif found_record
                    break
                elseif cs == 0
                    throw(EOFError())
                elseif p > p_eof ≥ 0
                    error("incomplete $($(reader_type)) input on line ", linenum)
                else
                    hits_eof = BufferedStreams.fillbuffer!(stream) == 0
                    p = stream.position
                    p_end = stream.available
                    if hits_eof
                        p_eof = p_end
                    end
                end
            end

            @assert isfilled(record)
            return record
        end
    end
end

end
