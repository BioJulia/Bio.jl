# BAM Overlap
# ===========

immutable OverlapIterator{T}
    reader::Reader{T}
    refname::String
    interval::UnitRange{Int}
end

function Base.iteratorsize{T}(::Type{OverlapIterator{T}})
    return Base.SizeUnknown()
end

function Base.eltype{T}(::Type{OverlapIterator{T}})
    return Record
end

function Bio.Intervals.eachoverlap(reader::Reader, interval::Bio.Intervals.Interval)
    return Bio.Intervals.eachoverlap(reader, interval.seqname, interval.first:interval.last)
end

function Bio.Intervals.eachoverlap(reader::Reader, refname::AbstractString, interval::UnitRange)
    return OverlapIterator(reader, String(refname), interval)
end


# Iterator
# --------

type OverlapIteratorState
    # reference index
    refindex::Int

    # possibly overlapping chunks
    chunks::Vector{Bio.Intervals.Chunk}

    # current chunk index
    chunkid::Int

    # pre-allocated record
    record::Record
end

function Base.start(iter::OverlapIterator)
    refindex = findfirst(iter.reader.refseqnames, iter.refname)
    if refindex == 0
        throw(ArgumentError("sequence name $(iter.refname) is not found in the header"))
    end
    @assert !isnull(iter.reader.index)
    chunks = Bio.Intervals.overlapchunks(get(iter.reader.index).index, refindex, iter.interval)
    if !isempty(chunks)
        seek(iter.reader, first(chunks).start)
    end
    return OverlapIteratorState(refindex, chunks, 1, Record())
end

function Base.done(iter::OverlapIterator, state)
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

function Base.next(::OverlapIterator, state)
    return copy(state.record), state
end

function isoverlapping(record::Record, refid_::Integer, interval::UnitRange)
    return refid(record) == refid_ && position(record) ≤ last(interval) && rightposition(record) ≥ first(interval)
end
