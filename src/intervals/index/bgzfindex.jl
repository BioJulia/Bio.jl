# BGZF Index
# ==========
#
# An index type for BGZFStream.
#
# The details of the internal is specified in
# https://samtools.github.io/hts-specs/SAMv1.pdf.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# binning index
const BinIndex = Dict{UInt32,Vector{Chunk}}

# linear index
const LinearIndex = Vector{VirtualOffset}

# Metadata providing a summary of the number of mappend/unmapped reads.
type PseudoBin
    # file range of unmapped reads
    unmapped::Chunk

    # number of mapped read segments
    n_mapped::Int64

    # number of unmapped read segments
    n_unmapped::Int64
end

# Index for BGZFStream; used in BAI and Tabix index.
type BGZFIndex
    # indexes of contigs (chromosomes)
    data::Vector{Tuple{BinIndex,LinearIndex,Nullable{PseudoBin}}}
end

# 16Kbp
const LinearWindowSize = 16 * 1024

# Find chunks overlapping with `(seqid, interval)` in `index`.
function overlapchunks(index::BGZFIndex, seqid::Integer, interval::UnitRange)
    if !(1 ≤ seqid ≤ endof(index.data))
        throw(ArgumentError("sequence id $(seqid) is out of range"))
    end

    if isempty(interval)
        return Chunk[]
    end

    binindex, linindex, pbin = index.data[seqid]
    bins = reg2bins(first(interval), last(interval))
    ret = Chunk[]
    idx = cld(first(interval), LinearWindowSize)
    if endof(linindex) ≥ idx
        # `linindex` may be empty for contigs with no records
        offset = linindex[idx]
        for bin in bins
            if haskey(binindex, bin)
                for chunk in binindex[bin]
                    if chunk.stop > offset
                        push!(ret, chunk)
                    end
                end
            end
        end
    end

    # tidy up the list of chunks
    sort!(ret)
    reduce!(ret)

    return ret
end

# Calculate bins overlapping a region [from, to] (one-based).
function reg2bins(from, to)
    bins = UInt32[]
    bin_start = 0
    for scale in 29:-3:14
        for k in ((from - 1) >> scale):((to - 1) >> scale)
            push!(bins, bin_start + k)
        end
        bin_start = 8 * bin_start + 1
    end
    return bins
end

# Merge chunks so as to minimize the number of seek operations.
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

# Read `n_refs` BAI/Tabix-compatible indexes from `input`.
function read_bgzfindex(input, n_refs)
    indexes = Tuple{BinIndex,LinearIndex,Nullable{PseudoBin}}[]
    for _ in 1:n_refs
        # load a binning index (and a pseudo bin)
        n_bins = read(input, Int32)
        binindex = BinIndex()
        pbin = Nullable{PseudoBin}()
        for _ in 1:n_bins
            bin = read(input, UInt32)
            n_chunks = read(input, Int32)
            if bin == 37450
                # pseudo-bin
                @assert n_chunks == 2
                chunk_beg = read(input, UInt64)
                chunk_end = read(input, UInt64)
                n_mapped = read(input, UInt64)
                n_unmapped = read(input, UInt64)
                pbin = Nullable(PseudoBin(
                    Chunk(chunk_beg, chunk_end),
                    n_mapped,
                    n_unmapped))
            else
                chunks = Chunk[]
                for i in 1:n_chunks
                    chunk_beg = read(input, UInt64)
                    chunk_end = read(input, UInt64)
                    push!(chunks, Chunk(chunk_beg, chunk_end))
                end
                binindex[bin] = chunks
            end
        end

        # load a linear index
        n_intvs = read(input, Int32)
        linindex = LinearIndex()
        for _ in 1:n_intvs
            push!(linindex, read(input, UInt64))
        end

        push!(indexes, (binindex, linindex, pbin))
    end
    return BGZFIndex(indexes)
end
