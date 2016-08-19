# BAI
# ===

# Binning index
typealias BinIndex Dict{UInt32,Vector{Chunk}}

# Linear index
typealias LinearIndex Vector{VirtualOffset}

const LinearWindowSize = 16 * 1024

type PseudoBin
    unmapped::Chunk
    n_mapped::UInt64
    n_unmapped::UInt64
end

type BAI
    indexes::Vector{Tuple{BinIndex,LinearIndex,Nullable{PseudoBin}}}
    n_no_coors::Nullable{UInt64}
end

function Base.read(input::IO, ::Type{BAI})
    # magic bytes
    B = read(input, UInt8)
    A = read(input, UInt8)
    I = read(input, UInt8)
    x = read(input, UInt8)
    if B != UInt8('B') || A != UInt8('A') || I != UInt8('I') || x != 0x01
        error("input is not a valid BAI file")
    end

    n_refs = read(input, Int32)
    indexes = read_indexes(input, n_refs)

    if !eof(input)
        n_no_coors = Nullable(read(input, UInt64))
    else
        n_no_coors = Nullable{UInt64}()
    end

    return BAI(indexes, n_no_coors)
end

# Read indexes for BAI and Tabix file formats.
function read_indexes(input, n_refs)
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

    return indexes
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
