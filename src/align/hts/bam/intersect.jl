# Intersect
# =========

type BAMIntersectionIterator
    reader::BAMReader
    chunks::Vector{Chunk}
    refindex::Int
    interval::UnitRange{Int}
end

Base.iteratorsize(::BAMIntersectionIterator) = Base.SizeUnknown()

function Base.start(iter::BAMIntersectionIterator)
    if !isempty(iter.chunks)
        seek(iter.reader, first(iter.chunks).start)
    end
    rec = BAMRecord()
    return advance!(iter, rec, 1)
end

function Base.done(iter::BAMIntersectionIterator, i_rec)
    i, _ = i_rec
    return i > endof(iter.chunks)
end

function Base.next(iter::BAMIntersectionIterator, i_rec)
    i, rec = i_rec
    ret = copy(rec)
    return ret, advance!(iter, rec, i)
end

function advance!(iter, rec, i)
    while i ≤ endof(iter.chunks)
        chunk = iter.chunks[i]
        while virtualoffset(iter.reader.stream) < chunk.stop
            read!(iter.reader, rec)
            if isoverlapping(rec, iter.refindex, iter.interval)
                return i, rec
            end
        end
        i += 1
        if i ≤ endof(iter.chunks)
            seek(iter.reader, iter.chunks[i].start)
        end
    end
    return i, rec
end

function Bio.Intervals.isoverlapping(rec, refindex_, interval)
    return ismapped(rec) &&
        refindex(rec) == refindex_ &&
        position(rec) ≤ last(interval) &&
        rightmost_position(rec) ≥ first(interval)
end

function Base.intersect(reader::BAMReader, interval::Interval)
    return intersect(reader, interval.seqname, interval.first:interval.last)
end

function Base.intersect(reader::BAMReader, refname::AbstractString, interval::UnitRange)
    i = findfirst(reader.refseqnames, refname)
    if i == 0
        error("sequence name $refname is not in the header")
    end
    return intersect(reader, i, interval)
end

function Base.intersect(reader::BAMReader, refindex::Integer, interval::UnitRange)
    if isnull(reader.index)
        error("no index")
    end
    chunks = overlapchunks(get(reader.index).index, refindex, interval)
    return BAMIntersectionIterator(reader, chunks, refindex, interval)
end
