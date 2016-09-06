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
        return new(Ragel.State(fastaparser_start, input),
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
    reader.state.cs = fastaparser_start
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
