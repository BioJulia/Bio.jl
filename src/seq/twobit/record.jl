# 2bit Record
# ===========

type Record
    filled::Bool
    dnasize::UInt32
    blockcount::UInt32
    blockstarts::Vector{UInt32}
    blocksizes::Vector{UInt32}
    maskedblockcount::UInt32
    maskedblockstarts::Vector{UInt32}
    maskedblocksizes::Vector{UInt32}
    reserved::UInt32
    packeddna::Vector{UInt8}
end

"""
    TwoBit.Record()

Create an unfilled 2bit record.
"""
function Record()
    return Record(
        false,
        # dnasize-blocksizes
        0, 0, UInt32[], UInt32[],
        # maskedblockcount-reserved
        0, UInt32[], UInt32[], 0,
        # packeddna
        UInt8[])
end

function initialize!(record::Record)
    record.filled = false
    record.dnasize = 0
    record.blockcount = 0
    empty!(record.blockstarts)
    empty!(record.blocksizes)
    record.maskedblockcount = 0
    empty!(record.maskedblockstarts)
    empty!(record.maskedblocksizes)
    record.reserved = 0
    return record
end

function isfilled(record::Record)
    return record.filled
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        return (
            record1.dnasize           == record2.dnasize           &&
            record1.blockcount        == record2.blockcount        &&
            record1.blockstarts       == record2.blockstarts       &&
            record1.blocksizes        == record2.blocksizes        &&
            record1.maskedblockcount  == record2.maskedblockcount  &&
            record1.maskedblockstarts == record2.maskedblockstarts &&
            record1.maskedblocksizes  == record2.maskedblocksizes  &&
            record1.reserved          == record2.reserved          &&
            memcmp(pointer(record1.packeddna), pointer(record2.packeddna), cld(record1.dnasize, 4)) == 0)
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "         length: ", record.dnasize)
        println(io, "       sequence: ", truncate(sequence(record), 40))
        println(io, "       N blocks: ", record.blockcount)
          print(io, "  masked blocks: ", record.maskedblockcount)
    else
        print(io, " <not filled>")
    end
end

function truncate(seq, len)
    if length(seq) > len
        return "$(seq[1:len-1])…"
    else
        return string(seq)
    end
end


# Accessors
# ---------

function hassequence(record::Record)
    return isfilled(record)
end

"""
    sequence([::Type{S},] record::Record)::S

Get the sequence of `record` as `S`.

`S` can be either `Bio.Seq.ReferenceSequence` or `Bio.Seq.DNASequence`.
If `S` is omitted, the default type is `Bio.Seq.ReferenceSequence`.
"""
function sequence(record::Record)
    return sequence(Bio.Seq.ReferenceSequence, record)
end

function sequence(::Type{Bio.Seq.ReferenceSequence}, record::Record)
    checkfilled(record)
    data = decode_sequence(record.packeddna, record.dnasize, 2, twobit2refseq_table)
    nmask = falses(record.dnasize)
    for i in 1:record.blockcount
        nmask[record.blockstarts[i] + (1:record.blocksizes[i])] = true
    end
    return Bio.Seq.ReferenceSequence(data, nmask, 1:record.dnasize)
end

function sequence(::Type{Bio.Seq.DNASequence}, record::Record)
    checkfilled(record)
    data = decode_sequence(record.packeddna, record.dnasize, 4, twobit2dnaseq_table)
    seq = Bio.Seq.DNASequence(data, 1:record.dnasize, false)
    for i in 1:record.blockcount
        seq[record.blockstarts[i] + (1:record.blocksizes[i])] = Bio.Seq.DNA_N
    end
    return seq
end

function Bio.sequence(record::Record)
    return sequence(record)
end

function Bio.hassequence(record::Record)
    return hassequence(record)
end

"""
    maskedblocks(record::Record)::Vector{UnitRange{Int}}

Get the masked blocks.
"""
function maskedblocks(record::Record)
    checkfilled(record)
    blocks = Vector{UnitRange{Int}}(record.maskedblockcount)
    for i in 1:endof(blocks)
        blocks[i] = record.maskedblockstarts[i] + (1:record.maskedblocksizes[i])
    end
    return blocks
end

function decode_sequence(packeddna, seqlen, nbits, table)
    @assert nbits ∈ (2, 4)
    data = zeros(UInt64, cld(seqlen, div(64, nbits)))
    stop = Bio.Seq.BitIndex(seqlen, nbits)
    i = Bio.Seq.BitIndex(1, nbits)
    j = 1
    while i ≤ stop
        data[Bio.Seq.index(i)] |= table[Int(packeddna[j])+1] << Bio.Seq.offset(i)
        i += 4 * nbits
        j += 1
    end
    return data
end

const twobit2refseq_table = let
    # T: 00, C: 01, A: 10, G: 11
    f(x) = x == 0b00 ? UInt64(3) :
           x == 0b01 ? UInt64(1) :
           x == 0b10 ? UInt64(0) :
           x == 0b11 ? UInt64(2) : error()
    tcag = 0b00:0b11
    tbl = UInt64[]
    for x in tcag, y in tcag, z in tcag, w in tcag
        push!(tbl, f(x) | f(y) << 2 | f(z) << 4 | f(w) << 6)
    end
    tbl
end

const twobit2dnaseq_table = let
    # T: 00, C: 01, A: 10, G: 11
    f(x) = x == 0b00 ? UInt64(Bio.Seq.DNA_T) :
           x == 0b01 ? UInt64(Bio.Seq.DNA_C) :
           x == 0b10 ? UInt64(Bio.Seq.DNA_A) :
           x == 0b11 ? UInt64(Bio.Seq.DNA_G) : error()
    tcag = 0b00:0b11
    tbl = UInt64[]
    for x in tcag, y in tcag, z in tcag, w in tcag
        push!(tbl, f(x) | f(y) << 4 | f(z) << 8 | f(w) << 12)
    end
    tbl
end

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled 2bit record"))
    end
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end
