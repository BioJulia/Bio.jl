# FASTA Reader
# ============

"""
    FASTAReader(input::IO; index=nothing)
    FASTAReader{S}(input::IO; index=nothing)

Create a data reader of the FASTA file format.

When type parameter `S` is specified, the reader reads sequences in that type;
otherwise the reader tries to infer the sequence type based on the frequencies
of characters from the input.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *fai* is supported)
"""
type FASTAReader{S<:Sequence} <: Bio.IO.AbstractReader
    state::Ragel.State
    seqbuf::BufferedOutputStream{BufferedStreams.EmptyStream}
    index::Nullable{FASTAIndex}

    function FASTAReader(input::BufferedInputStream, index)
        return new(Ragel.State(fasta_machine.start_state, input),
                   BufferedOutputStream(), index)
    end
end

function Bio.IO.stream(reader::FASTAReader)
    return reader.state.stream
end

function Base.eltype{S}(::Type{FASTAReader{S}})
    return FASTASeqRecord{S}
end

function FASTAReader(input::IO; index=nothing)
    # alphabet type will be inferred from data
    return FASTAReader{BioSequence}(input, index=index)
end

function (::Type{FASTAReader{S}}){S<:Sequence}(input::IO; index=nothing)
    if isa(index, AbstractString)
        index = FASTAIndex(index)
    else
        if index != nothing
            error("unrecognizable index argument")
        end
    end
    return FASTAReader{S}(BufferedInputStream(input), index)
end

function Base.getindex(reader::FASTAReader, name::AbstractString)
    if isnull(reader.index)
        error("no index")
    end
    seekrecord(reader.state.stream, get(reader.index), name)
    reader.state.cs = fasta_machine.start_state
    return read(reader)
end

# Predict sequence type based on character frequencies in `seq[start:stop]`.
function predict(seq::Vector{UInt8}, start, stop)
    # count characters
    a = c = g = t = u = n = alpha = 0
    for i in start:stop
        @inbounds x = seq[i]
        if x == 0x41 || x == 0x61
            a += 1
        elseif x == 0x43 || x == 0x63
            c += 1
        elseif x == 0x47 || x == 0x67
            g += 1
        elseif x == 0x54 || x == 0x74
            t += 1
        elseif x == 0x55 || x == 0x75
            u += 1
        elseif x == 0x4e || x == 0x6e
            n += 1
        end
        if 0x41 ≤ x ≤ 0x5a || 0x61 ≤ x ≤ 0x7a
            alpha += 1
            if alpha ≥ 300 && t + u > 0 && a + c + g + t + u + n == alpha
                # pretty sure that the sequence is either DNA or RNA
                break
            end
        end
    end

    # the threshold (= 0.95) is somewhat arbitrary
    if (a + c + g + t + u + n) / alpha > 0.95
        if t ≥ u
            return DNASequence
        else
            return RNASequence
        end
    else
        return AminoAcidSequence
    end
end

info("compiling old FASTA")
const fasta_machine = (function ()
    re = Automa.RegExp

    lf          = re"\n"
    newline     = re"\r?" * lf
    hspace      = re"[ \t\v]"
    whitespace  = re.space() | newline
    identifier  = re.rep1(re.any() \ re.space())
    description = (re.any() \ hspace) * re"[^\r\n]*"
    letters     = re.rep1(re.any() \ re.space() \ re">")
    sequence    = re.cat(re.rep(whitespace), re.opt(letters), re.rep(re.rep1(whitespace) * letters))
    record      = re.cat('>', identifier, re.opt(re.rep1(hspace) * description), newline, sequence, re.rep(whitespace))
    fasta       = re.rep(whitespace) * re.rep(record)

    lf.actions[:enter]          = [:count_line]
    identifier.actions[:enter]  = [:mark]
    identifier.actions[:exit]   = [:identifier]
    description.actions[:enter] = [:mark]
    description.actions[:exit]  = [:description]
    letters.actions[:enter]     = [:mark]
    letters.actions[:exit]      = [:letters]
    record.actions[:exit]       = [:record]

    return Automa.compile(fasta)
end)()

@inline function anchor!(stream, p)
    stream.anchor = p
end

@inline function upanchor!(stream)
    @assert stream.anchor != 0 "upanchor! called with no anchor set"
    anchor = stream.anchor
    stream.anchor = 0
    return anchor
end

const fasta_actions = Dict(
    :count_line  => :(linenum += 1),
    :mark        => :(anchor!(stream, p)),
    :identifier  => :(copy!(output.name, stream.buffer, upanchor!(stream), p - 1)),
    :description => :(copy!(output.metadata.description, stream.buffer, upanchor!(stream), p - 1)),
    :letters     => :(append!(reader.seqbuf, stream.buffer, upanchor!(stream), p - 1)),
    :record      => :(found_record = true; @escape))

function Base.read!(reader::FASTAReader, output::FASTASeqRecord)
    return _read!(reader, reader.state, output)
end

@eval function _read!(reader::FASTAReader, state::Ragel.State, output::FASTASeqRecord)
    cs = state.cs
    linenum = state.linenum
    stream = state.stream
    data = stream.buffer
    p = stream.position
    p_end = stream.available
    p_eof = -1
    found_record = false

    while true
        $(Automa.generate_exec_code(fasta_machine, actions=fasta_actions, code=:goto, check=false))

        state.cs = cs
        state.finished = cs == 0
        state.linenum = linenum
        stream.position = p

        if cs < 0
            error("FASTA file format error on line ", linenum)
        elseif found_record
            if seqtype(typeof(output)) == ReferenceSequence
                output.seq = ReferenceSequence(reader.seqbuf.buffer, 1, length(reader.seqbuf))
            elseif seqtype(typeof(output)) == BioSequence
                ET = predict(reader.seqbuf.buffer, 1, length(reader.seqbuf))
                if ET == typeof(output.seq)
                    resize!(output.seq, length(reader.seqbuf))
                    encode_copy!(output.seq, 1, reader.seqbuf.buffer, 1, length(reader.seqbuf))
                else
                    output.seq = ET(reader.seqbuf.buffer, 1, length(reader.seqbuf))
                end
            else
                resize!(output.seq, length(reader.seqbuf))
                encode_copy!(output.seq, 1, reader.seqbuf.buffer, 1, length(reader.seqbuf))
            end
            empty!(reader.seqbuf)
            break
        elseif cs == 0
            throw(EOFError())
        elseif p > p_eof ≥ 0
            error("incomplete FASTA input on line ", linenum)
        else
            hits_eof = BufferedStreams.fillbuffer!(stream) == 0
            p = stream.position
            p_end = stream.available
            if hits_eof
                p_eof = p_end
            end
        end
    end

    return output
end
