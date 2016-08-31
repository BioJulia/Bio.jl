# 2bit format
# ===========
#
# Reader and writer of the 2bit file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
The 2bit file format stores multiple DNA sequences (up to 4 Gbp total) as a
compact randomly-accessible format.
See https://genome.ucsc.edu/FAQ/FAQformat.html#format7 for the details.
"""
immutable TwoBit <: Bio.IO.FileFormat end

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

function Base.open(filename::AbstractString, ::Type{TwoBit})
    return TwoBitReader(open(filename))
end

function Base.eltype(::Type{TwoBitReader})
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


# Reader
# ------

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


# Writer
# ------

type TwoBitWriter{T<:IO} <: Bio.IO.AbstractWriter
    # output stream
    output::T

    # sequence names
    names::Vector{String}

    # bit vector to check if each sequence is already written or not
    written::BitVector
end

function TwoBitWriter(output::IO, names::AbstractVector)
    writer = TwoBitWriter(output, names, falses(length(names)))
    write_header(writer)
    write_index(writer)
    return writer
end

function Bio.IO.stream(writer::TwoBitWriter)
    return writer.output
end

function Base.close(writer::TwoBitWriter)
    if !all(writer.written)
        error("one or more sequences are not written")
    end
    close(writer.output)
end

function write_header(writer::TwoBitWriter)
    output = writer.output
    n = 0
    n += write(output, 0x1A412743)
    n += write(output, UInt32(0))
    n += write(output, UInt32(length(writer.names)))
    n += write(output, UInt32(0))
    return n
end

function write_index(writer::TwoBitWriter)
    output = writer.output
    n = 0
    for name in writer.names
        n += write(output, UInt8(length(name)))
        n += write(output, name)
        n += write(output, UInt32(0))
    end
    return n
end

# Update the file offset of a sequence in the index section.
function update_offset(writer::TwoBitWriter, seqname, seqoffset)
    @assert seqname ∈ writer.names
    output = writer.output
    offset = 16
    for name in writer.names
        offset += sizeof(UInt8) + length(name)
        if name == seqname
            old = position(output)
            seek(output, offset)
            write(output, UInt32(seqoffset))
            seek(output, old)
            return
        end
        offset += sizeof(UInt32)
    end
end

function Base.write(writer::TwoBitWriter, record::SeqRecord)
    i = findfirst(writer.names, record.name)
    if i == 0
        error("sequence \"", record.name, "\" doesn't exist in the writing list")
    elseif writer.written[i]
        error("sequence \"", record.name, "\" is already written")
    end

    output = writer.output
    update_offset(writer, record.name, position(output))

    n = 0
    n += write(output, UInt32(length(record.seq)))
    n += write_n_blocks(output, record.seq)
    n += write_masked_blocks(output, record.metadata)
    n += write(output, UInt32(0))  # reserved bytes
    n += write_twobit_sequence(output, record.seq)

    writer.written[i] = true
    return n
end

function make_n_blocks(seq)
    starts = UInt32[]
    sizes = UInt32[]
    i = 1
    while i ≤ endof(seq)
        nt = seq[i]
        if nt == DNA_N
            start = i - 1  # 0-based index
            push!(starts, start)
            while i ≤ endof(seq) && seq[i] == DNA_N
                i += 1
            end
            push!(sizes, (i - 1) - start)
        elseif isambiguous(nt)
            error("ambiguous nucleotide except N is not supported")
        else
            i += 1
        end
    end
    return starts, sizes
end

function write_n_blocks(output, seq)
    blockstarts, blocksizes = make_n_blocks(seq)
    @assert length(blockstarts) == length(blocksizes)
    n = 0
    n += write(output, UInt32(length(blockstarts)))
    n += write(output, blockstarts)
    n += write(output, blocksizes)
    return n
end

function write_masked_blocks(output, metadata)
    n = 0
    if isa(metadata, Vector{UnitRange{Int}})
        n += write(output, UInt32(length(metadata)))
        for mblock in metadata
            n += write(output, UInt32(first(mblock) - 1))  # 0-based
        end
        for mblock in metadata
            n += write(output, UInt32(length(mblock)))
        end
    elseif metadata === nothing
        n += write(output, UInt32(0))
    else
        error("metadata is not serializable in the 2bit file format")
    end
    return n
end

function write_twobit_sequence(output, seq)
    n = 0
    i = 4
    while i ≤ endof(seq)
        x::UInt8 = 0
        x |= nuc2twobit(seq[i-3]) << 6
        x |= nuc2twobit(seq[i-2]) << 4
        x |= nuc2twobit(seq[i-1]) << 2
        x |= nuc2twobit(seq[i-0]) << 0
        n += write(output, x)
        i += 4
    end
    r = length(seq) % 4
    if r > 0
        let x::UInt8 = 0
            i = endof(seq) - r + 1
            while i ≤ endof(seq)
                x = x << 2 | nuc2twobit(seq[i])
                i += 1
            end
            x <<= (4 - r) * 2
            n += write(output, x)
        end
    end
    return n
end

function nuc2twobit(nt::DNANucleotide)
    return (
        nt == DNA_A ? 0b10 :
        nt == DNA_C ? 0b01 :
        nt == DNA_G ? 0b11 :
        nt == DNA_T ? 0b00 :
        nt == DNA_N ? 0b00 : error())
end
