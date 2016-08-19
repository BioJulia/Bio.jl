# Tabix
# =====
#
# Generic index for tab-delimited files.
#
# Li, Heng. "Tabix: fast retrieval of sequence features from generic
# TAB-delimited files." Bioinformatics 27.5 (2011): 718-719.
# Specification: http://samtools.github.io/hts-specs/tabix.pdf

include("index.jl")
include("reader.jl")

"""
    overlapchunks(index::Tabix, seqname::String, interval::UnitRange)
    overlapchunks(index::Tabix, seqid::Integer, interval::UnitRange)

Return chunks possibly overlapping the range specified by `seqname` (or `seqid`)
and `interval`; `seqid` and `interval` must be 1-based index and inclusive.

NOTE: Records within the returned chunks are not guaranteed to overlap the query
region. The user need to check the condition as:
    line = readline(stream)
    values = split(chomp(line), '\t')
    seqname = value[index.columns[1]]
    start = parse(Int, index.columns[2])
    stop = parse(Int, index.columns[3]) - 1  # needs -1 for the BED format
    if seqname == target_seqname && intersect(start:stop, target_interval)
        # the line overlaps the query region
    end
"""
function overlapchunks(index::Union{BAI,Tabix}, seqid::Integer, interval::UnitRange)
    if !(1 ≤ seqid ≤ endof(index.indexes))
        throw(ArgumentError("sequence id $(seqid) is out of range"))
    end

    if isempty(interval)
        return Chunk[]
    end

    binindex, linindex, pbin = index.indexes[seqid]
    bins = reg2bins(first(interval), last(interval))
    offset = linindex[cld(first(interval), LinearWindowSize)]
    ret = Chunk[]
    for bin in bins
        if haskey(binindex, bin)
            for chunk in binindex[bin]
                if chunk.stop > offset
                    push!(ret, chunk)
                end
            end
        end
    end

    sort!(ret)
    reduce!(ret)

    return ret
end

function overlapchunks(index::Tabix, seqname::AbstractString, interval::UnitRange)
    seqid = findfirst(index.names, seqname)
    if seqid == 0
        throw(ArgumentError("sequence name $(seqname) is not included in the index"))
    end
    return overlapchunks(index, seqid, interval)
end

# Merge chunks to minimize the number of seek operations.
function reduce!(chunks)
    @assert issorted(chunks)
    i = 1
    while i < endof(chunks)
        chunk = chunks[i]
        next = chunks[i+1]
        if chunk.stop < next.start
            # neither overlapping nor adjacent
            i += 1
            continue
        end
        chunks[i] = Chunk(chunk.start, max(chunk.stop, next.stop))
        deleteat!(chunks, i + 1)
    end
    return chunks
end
