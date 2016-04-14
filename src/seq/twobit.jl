# 2bit format
# ===========

"""
.2bit format

.2bit file format stores multiple DNA sequences (up to 4 Gib total) in a compact
randomly-accessible format: https://genome.ucsc.edu/FAQ/FAQformat.html#format7.
"""
immutable TwoBit <: FileFormat end

immutable TwoBitParser <: AbstractParser
    input::IO
    names::Vector{ASCIIString}
    offsets::Vector{UInt32}
    swap::Bool
end

function TwoBitParser(input::IO)
    names, offsets, swap = loadheader(input)
    return TwoBitParser(input, names, offsets, swap)
end

function Base.open(filename::AbstractString, ::Type{TwoBit})
    input = open(filename)
    finalizer(input, close)
    return TwoBitParser(input)
end

# iterator
Base.eltype(::Type{TwoBitParser}) = SeqRecord{DNASequence,Vector{UnitRange{Int}}}
Base.length(p::TwoBitParser) = length(p.names)
Base.start(p::TwoBitParser) = 1
Base.done(p::TwoBitParser, k::Int) = k > endof(p)
Base.next(p::TwoBitParser, k::Int) = p[k], k + 1


# Indexing
# --------

function Base.checkbounds(p::TwoBitParser, k::Integer)
    if 1 ≤ k ≤ endof(p)
        return true
    end
    throw(BoundsError(p, k))
end

function Base.checkbounds(p::TwoBitParser, r::Range)
    if isempty(r) || (1 ≤ first(r) && last(r) ≤ endof(p))
        return true
    end
    throw(BoundsError(p, k))
end

Base.endof(p::TwoBitParser) = length(p)

function Base.getindex(p::TwoBitParser, k::Integer)
    checkbounds(p, k)
    seek(p.input, p.offsets[k])
    seq, mblocks = readseq(p.input, p.swap)
    return SeqRecord{DNASequence,Vector{UnitRange{Int}}}(p.names[k], seq, mblocks)
end

function Base.getindex(p::TwoBitParser, name::AbstractString)
    k = findfirst(p.names, name)
    if k == 0
        throw(KeyError(name))
    end
    return p[k]
end

function Base.getindex(p::TwoBitParser, r::Range)
    checkbounds(p, r)
    return [p[k] for k in r]
end

function Base.getindex{S<:AbstractString}(p::TwoBitParser, names::AbstractVector{S})
    return [p[name] for name in names]
end


# Parser
# ------

function loadheader(io::IO)
    # header
    signature = read(io, UInt32)
    if signature == 0x1A412743
        swap = false
    elseif signature == 0x4327411A
        swap = true
    else
        error("invalid signature")
    end

    version  = read32(io, swap)
    seqcount = read32(io, swap)
    reserved = read32(io, swap)
    @assert version == reserved == 0

    # file index
    names = ASCIIString[]
    offsets = UInt32[]
    for i in 1:seqcount
        namesize = read(io, UInt8)
        name = read(io, UInt8, Int(namesize))
        offset = read32(io, swap)
        push!(names, name)
        push!(offsets, offset)
    end

    return names, offsets, swap
end

# mapping table from .2bit DNA encoding to Bio.jl DNA encoding
const twobit_decode_table = let
    # T: 00, C: 01, A: 10, G: 11
    f(x) = x == 0b00 ? UInt64(DNA_T) :
           x == 0b01 ? UInt64(DNA_C) :
           x == 0b10 ? UInt64(DNA_A) :
           x == 0b11 ? UInt64(DNA_G) : error()
    tcag = 0b00:0b11
    tbl = UInt64[]
    for a in tcag, b in tcag, c in tcag, d in tcag
        x = UInt64(0)
        x |= f(a)
        x |= f(b) <<  4
        x |= f(c) <<  8
        x |= f(d) << 12
        push!(tbl, x)
    end
    tbl
end

twobit_decode(pack::UInt8) = twobit_decode_table[UInt(pack)+1]

# read 32 bits and swap bytes if necessary
function read32(io, swap)
    x = read(io, UInt32)
    if swap
        x = bswap(x)
    end
    return x
end

function read32(io, n, swap)
    xs = read(io, UInt32, n)
    if swap
        map!(bswap, xs)
    end
    return xs
end

function readseq(io, swap)
    # read record header
    dnasize::Int = read32(io, swap)
    blockcount::Int = read32(io, swap)  # N blocks
    blockstarts = read32(io, blockcount, swap)
    blocksizes  = read32(io, blockcount, swap)
    mblockcount::Int = read32(io, swap)  # masked blocks
    mblockstarts = read32(io, mblockcount, swap)
    mblocksizes  = read32(io, mblockcount, swap)
    reserved = read32(io, swap)
    @assert reserved == 0

    # read packed DNA
    data = zeros(UInt64, cld(dnasize, 16))
    j = 1
    r = 0
    for _ in 1:cld(dnasize, 4)
        pack = read(io, UInt8)
        #data[j] |= twobit_decode_table[UInt(pack)+1] << r
        data[j] |= twobit_decode(pack) << r
        r += 16
        if r ≥ 64
            j += 1
            r = 0
        end
    end

    # fill N blocks
    seq = DNASequence(data, 1:dnasize, false)
    for k in 1:blockcount
        s = blockstarts[k]
        l = blocksizes[k]
        for i in s+1:s+l
            seq[i] = DNA_N
        end
    end

    # make masked blocks
    mblocks = UnitRange{Int}[]
    for k in 1:mblockcount
        s = mblockstarts[k]
        l = mblocksizes[k]
        push!(mblocks, s+1:s+l)
    end

    return seq, mblocks
end


# Serializer
# ----------

function writeformat(io::IO, ::Type{TwoBit}, records)
    # compute N blocks
    blockstarts = Vector{UInt32}[]
    blocksizes  = Vector{UInt32}[]
    for r in records
        seq = r.seq
        starts = UInt32[]
        sizes  = UInt32[]
        i = 1
        while i ≤ endof(seq)
            if seq[i] == DNA_N
                start = i - 1  # 0-based index
                push!(starts, start)
                while i ≤ endof(seq) && seq[i] == DNA_N
                    i += 1
                end
                push!(sizes, (i - 1) - start)
            else
                if seq[i] > DNA_T
                    error("ambiguous nucleotide except N is not supported")
                end
                i += 1
            end
        end
        push!(blockstarts, starts)
        push!(blocksizes, sizes)
    end
    @assert length(blockstarts) == length(blocksizes)

    # write header
    write(io, 0x1A412743)
    write(io, UInt32(0))
    write(io, UInt32(length(records)))
    write(io, UInt32(0))

    # write index
    offset = 0  # keep track of the file offset for each sequence record
    offset += 16  # header size
    offset += sum([1 + length(r.name) + 4 for r in records])  # index size
    for (k, r) in enumerate(records)
        write(io, UInt8(length(r.name)))
        write(io, r.name)
        write(io, UInt32(offset))

        # increment the offset by the size of the sequence record
        offset += 4                               # dna size
        offset += 4                               # block count
        offset += 4 * length(blockstarts[k]) * 2  # blockstarts and blocksizes
        offset += 4                               # masked block count
        offset += 4 * length(r.metadata) * 2      # masked blockstarts and blocksizes
        offset += 4                               # reserved bytes
        offset += cld(length(r.seq), 4)           # DNA sequence
    end

    # write sequence records
    for (k, r) in enumerate(records)
        seq = r.seq
        write(io, UInt32(length(seq)))

        # write N blocks
        write(io, UInt32(length(blockstarts[k])))
        write(io, blockstarts[k])
        write(io, blocksizes[k])

        # write masked blocks
        masks = r.metadata
        @assert eltype(masks) <: UnitRange
        write(io, UInt32(length(masks)))
        for m in masks
            write(io, UInt32(first(m) - 1))  # 0-based index
            write(io, UInt32(length(m)))
        end

        # write reserved bytes
        write(io, UInt32(0))

        # write packed DNA sequence
        i = 4
        while i ≤ endof(seq)
            x::UInt8 = 0
            x |= twobit_encode(seq[i-3]) << 6
            x |= twobit_encode(seq[i-2]) << 4
            x |= twobit_encode(seq[i-1]) << 2
            x |= twobit_encode(seq[i])
            write(io, x)
            i += 4
        end
        r = length(seq) % 4
        if r > 0
            y::UInt8 = 0
            i = endof(seq) - r + 1
            while i ≤ endof(seq)
                y = y << 2 | twobit_encode(seq[i])
                i += 1
            end
            y <<= (4 - r) * 2
            write(io, y)
        end
    end
end

function twobit_encode(nt::DNANucleotide)
    return (
        nt == DNA_A ? 0b10 :
        nt == DNA_C ? 0b01 :
        nt == DNA_G ? 0b11 :
        nt == DNA_T ? 0b00 :
        nt == DNA_N ? 0b00 : error()
    )
end
