# 2bit Reader
# ===========

"""
    TwoBitReader(input::IO)

Create a data reader of the 2bit file format.

# Arguments
* `input`: data source
"""
type TwoBitReader{T<:IO} <: Bio.IO.AbstractReader
    # input stream
    input::T

    # sequence names
    names::Vector{String}

    # file offsets
    offsets::Vector{UInt32}

    # byte-swapped or not
    swapped::Bool
end

function TwoBitReader(input::IO)
    swapped, seqcount = read_2bit_header(input)
    names, offsets = read_2bit_index(input, swapped, seqcount)
    @assert seqcount == length(names) == length(offsets)
    return TwoBitReader(input, names, offsets, swapped)
end

function Bio.IO.stream(reader::TwoBitReader)
    return reader.input
end

function Base.eltype{S}(::Type{TwoBitReader{S}})
    return SeqRecord{ReferenceSequence,Vector{UnitRange{Int}}}
end

function Base.length(p::TwoBitReader)
    return length(p.names)
end

function Base.start(p::TwoBitReader)
    return 1
end

function Base.done(p::TwoBitReader, k)
    return k > endof(p)
end

function Base.next(p::TwoBitReader, k)
    return p[k], k + 1
end


# Random access
# -------------

function Base.checkbounds(p::TwoBitReader, k::Integer)
    if 1 ≤ k ≤ endof(p)
        return true
    end
    throw(BoundsError(p, k))
end

function Base.checkbounds(p::TwoBitReader, r::Range)
    if isempty(r) || (1 ≤ first(r) && last(r) ≤ endof(p))
        return true
    end
    throw(BoundsError(p, k))
end

function Base.endof(p::TwoBitReader)
    return length(p)
end

function Base.getindex(p::TwoBitReader, k::Integer)
    checkbounds(p, k)
    seek(p.input, p.offsets[k])
    seq, mblocks = read_2bit_seq(p.input, p.swapped)
    return SeqRecord{ReferenceSequence,Vector{UnitRange{Int}}}(p.names[k], seq, mblocks)
end

function Base.getindex(p::TwoBitReader, name::AbstractString)
    k = findfirst(p.names, name)
    if k == 0
        throw(KeyError(name))
    end
    return p[k]
end

function Base.getindex(p::TwoBitReader, r::Range)
    checkbounds(p, r)
    return [p[k] for k in r]
end

function Base.getindex{S<:AbstractString}(p::TwoBitReader, names::AbstractVector{S})
    return [p[name] for name in names]
end

function read_2bit_header(input)
    signature = read(input, UInt32)
    if signature == 0x1A412743
        swapped = false
    elseif signature == 0x4327411A
        swapped = true
    else
        error("invalid 2bit signature")
    end
    version  = read32(input, swapped)
    seqcount = read32(input, swapped)
    reserved = read32(input, swapped)
    @assert version == reserved == 0
    return swapped, seqcount
end

function read_2bit_index(input, swapped, seqcount)
    names = String[]
    offsets = UInt32[]
    for i in 1:seqcount
        namesize = read(input, UInt8)
        name = read(input, UInt8, namesize)
        offset = read32(input, swapped)
        push!(names, String(name))
        push!(offsets, offset)
    end
    return names, offsets
end

function read_2bit_seq(input, swap)
    # read record header
    dnasize = read32(input, swap)
    blockcount = read32(input, swap)  # N blocks
    blockstarts = read32(input, blockcount, swap)
    blocksizes = read32(input, blockcount, swap)
    mblockcount = read32(input, swap)  # masked blocks
    mblockstarts = read32(input, mblockcount, swap)
    mblocksizes = read32(input, mblockcount, swap)
    reserved = read32(input, swap)
    @assert reserved == 0

    # read packed DNAs
    data = zeros(UInt64, cld(dnasize, 32))
    i_stop = BitIndex(dnasize, 2)
    i = BitIndex(1, 2)
    while i ≤ i_stop
        x = read(input, UInt8)
        data[index(i)] |= twobit2refseq(x) << offset(i)
        i += 8
    end

    # make an N vector
    nmask = falses(Int(dnasize))
    for j in 1:blockcount
        s = blockstarts[j]
        for k in s+1:s+blocksizes[j]
            nmask[k] = true
        end
    end

    # make masked blocks
    mblocks = UnitRange{Int}[]
    for j in 1:mblockcount
        s = mblockstarts[j]
        push!(mblocks, s+1:s+mblocksizes[j])
    end

    return ReferenceSequence(data, nmask, 1:dnasize), mblocks
end

# mapping table from .2bit DNA encoding to Bio.jl DNA encoding
const twobit2refseq_table = let
    # T: 00, C: 01, A: 10, G: 11
    f(x) = x == 0b00 ? UInt64(3) :
           x == 0b01 ? UInt64(1) :
           x == 0b10 ? UInt64(0) :
           x == 0b11 ? UInt64(2) : error()
    tcag = 0b00:0b11
    tbl = UInt64[]
    for a in tcag, b in tcag, c in tcag, d in tcag
        x::UInt64 = 0
        x |= f(a) << 0
        x |= f(b) << 2
        x |= f(c) << 4
        x |= f(d) << 6
        push!(tbl, x)
    end
    tbl
end

twobit2refseq(x::UInt8) = twobit2refseq_table[Int(x)+1]

# read 32 bits and swap bytes if necessary
function read32(input, swap)
    x = read(input, UInt32)
    if swap
        x = bswap(x)
    end
    return x
end

function read32(input, n, swap)
    xs = read(input, UInt32, n)
    if swap
        map!(bswap, xs)
    end
    return xs
end
