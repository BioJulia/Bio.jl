# Intersect
# =========

immutable BAMIntersectionIterator{T}
    # data reader
    reader::BAMReader{T}

    # query interval
    refname::AbstractString
    interval::UnitRange{Int}
end

function Base.iteratorsize(::BAMIntersectionIterator)
    return Base.SizeUnknown()
end

function Base.eltype{T}(::Type{BAMIntersectionIterator{T}})
    return BAMRecord
end

function Base.show(io::IO, iter::BAMIntersectionIterator)
    print(
        io,
        summary(iter),
        "(<$(iter.refname):$(first(iter.interval))-$(last(iter.interval))>)")
end

type BAMIntersectionIteratorState
    # reference index
    refindex::Int

    # possibly overlapping chunks
    chunks::Vector{Bio.Intervals.Chunk}

    # current chunk index
    chunkid::Int

    # pre-allocated record
    record::BAMRecord
end

function Base.start(iter::BAMIntersectionIterator)
    refindex = findfirst(iter.reader.refseqnames, iter.refname)
    if refindex == 0
        throw(ArgumentError("sequence name $(iter.refname) is not found in the header"))
    end
    @assert !isnull(iter.reader.index)
    chunks = Bio.Intervals.overlapchunks(get(iter.reader.index).index, refindex, iter.interval)
    if !isempty(chunks)
        seek(iter.reader, first(chunks).start)
    end
    return BAMIntersectionIteratorState(refindex, chunks, 1, BAMRecord())
end

function Base.done(iter::BAMIntersectionIterator, state)
    while state.chunkid ≤ endof(state.chunks)
        chunk = state.chunks[state.chunkid]
        while BGZFStreams.virtualoffset(iter.reader.stream) < chunk.stop
            read!(iter.reader, state.record)
            if isoverlapping(state.record, state.refindex, iter.interval)
                return false
            end
        end
        state.chunkid += 1
        if state.chunkid ≤ endof(state.chunks)
            seek(iter.reader, state.chunks[state.chunkid].start)
        end
    end
    return true
end

function Base.next(iter::BAMIntersectionIterator, state)
    return copy(state.record), state
end

function Bio.isoverlapping(
        rec::BAMRecord,
        refindex_::Integer,
        interval::UnitRange)
    return refindex(rec) == refindex_ &&
        leftposition(rec) ≤ last(interval) &&
        rightposition(rec) ≥ first(interval)
end

function Base.intersect(reader::BAMReader, interval::Bio.Intervals.Interval)
    return intersect(reader, interval.seqname, interval.first:interval.last)
end

function Base.intersect(
        reader::BAMReader,
        refname::AbstractString,
        interval::UnitRange)
    return BAMIntersectionIterator(reader, refname, interval)
end
