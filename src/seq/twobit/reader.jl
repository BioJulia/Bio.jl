# 2bit Reader
# ===========

type Reader{T<:IO} <: Bio.IO.AbstractReader
    # input stream
    input::T
    # sequence names
    names::Vector{String}
    # file offsets
    offsets::Vector{UInt32}
    # byte-swapped or not
    swap::Bool
end

"""
    TwoBit.Reader(input::IO)

Create a data reader of the 2bit file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO)
    seqcount, swap = readheader!(input)
    names, offsets = readindex!(input, seqcount, swap)
    @assert seqcount == length(names) == length(offsets)
    return Reader(input, names, offsets, swap)
end

function Bio.IO.stream(reader::Reader)
    return reader.input
end

function Base.eltype{T}(::Type{Reader{T}})
    return Record
end

function Base.length(reader::Reader)
    return length(reader.names)
end

function Base.start(reader::Reader)
    return 1
end

function Base.done(reader::Reader, i)
    return i > endof(reader.names)
end

function Base.next(reader::Reader, i)
    return reader[i], i + 1
end

"""
    seqnames(reader::Reader)::Vector{String}

Get the sequence names.

Sequences are stored in this order.
"""
function seqnames(reader::Reader)
    return copy(reader.names)
end


# Random access
# -------------

function Base.checkbounds(reader::Reader, i::Integer)
    if 1 ≤ i ≤ endof(reader)
        return true
    end
    throw(BoundsError(reader, i))
end

function Base.checkbounds(reader::Reader, r::Range)
    if isempty(r) || (1 ≤ first(r) && last(r) ≤ endof(reader))
        return true
    end
    throw(BoundsError(reader, r))
end

function Base.endof(reader::Reader)
    return length(reader)
end

function Base.getindex(reader::Reader, i::Integer)
    checkbounds(reader, i)
    seek(reader.input, reader.offsets[i])
    record = Record()
    read!(reader, record)
    return record
end

function Base.getindex(reader::Reader, r::Range)
    checkbounds(reader, r)
    return [reader[i] for i in r]
end

function Base.getindex(reader::Reader, name::AbstractString)
    i = findfirst(reader.names, name)
    if i == 0
        throw(KeyError(name))
    end
    return reader[i]
end

function Base.getindex{S<:AbstractString}(reader::Reader, names::AbstractVector{S})
    return [reader[name] for name in names]
end


# Readers
# -------

const SIGNATURE = UInt32(0x1A412743)

function readheader!(input)
    signature = read(input, UInt32)
    if signature == SIGNATURE
        swap = false
    elseif signature == bswap(SIGNATURE)
        swap = true
    else
        error("invalid 2bit signature")
    end
    version  = read32(input, swap)
    seqcount = read32(input, swap)
    reserved = read32(input, swap)
    @assert version == reserved == 0
    return seqcount, swap
end

function readindex!(input, seqcount, swap)
    names = String[]
    offsets = UInt32[]
    for i in 1:seqcount
        namesize = read(input, UInt8)
        name = read(input, UInt8, namesize)
        offset = read32(input, swap)
        push!(names, String(name))
        push!(offsets, offset)
    end
    return names, offsets
end

function Base.read!(reader::Reader, record::Record)
    initialize!(record)
    record.dnasize = read32(reader)
    record.blockcount = read32(reader)
    read32!(reader, record.blockstarts, record.blockcount)
    read32!(reader, record.blocksizes, record.blockcount)
    record.maskedblockcount = read32(reader)
    read32!(reader, record.maskedblockstarts, record.maskedblockcount)
    read32!(reader, record.maskedblocksizes, record.maskedblockcount)
    record.reserved = read32(reader)
    datalen = cld(record.dnasize, 4)
    if length(record.packeddna) < datalen
        resize!(record.packeddna, datalen)
    end
    unsafe_read(reader.input, pointer(record.packeddna), datalen)
    record.filled = true
    return record
end

# read 32 bits and swap bytes if necessary
function read32(reader::Reader)
    return read32(reader.input, reader.swap)
end

function read32!(reader::Reader, out, n)
    resize!(out, n)
    read!(reader.input, out)
    if reader.swap
        map!(bswap, out)
    end
    return out
end

function read32(input::IO, swap)
    x = read(input, UInt32)
    if swap
        x = bswap(x)
    end
    return x
end
