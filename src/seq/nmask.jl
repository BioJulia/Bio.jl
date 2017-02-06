# NMask
# =====
#
# Compressed bit vector of N nucleotides.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable NMask
    blockmask::IndexableBitVectors.SucVector
    blocks::Vector{UInt64}
    len::Int
end

function NMask()
    return NMask(IndexableBitVectors.SucVector(), UInt64[], 0)
end

Base.copy(nmask::NMask) = NMask(copy(nmask.blockmask), copy(nmask.blocks), nmask.len)

function Base.convert(::Type{NMask}, bv::BitVector)
    n = length(bv.chunks)
    blockmask = BitVector()
    blocks = Vector{UInt64}()
    for i in 1:n
        chunk = bv.chunks[i]
        if chunk == 0
            # no N in this block
            push!(blockmask, false)
        else
            push!(blockmask, true)
            push!(blocks, chunk)
        end
    end
    return NMask(blockmask, blocks, length(bv))
end

function Base.convert(::Type{BitVector}, nmask::NMask)
    bv = BitVector()
    for i in 1:length(nmask)
        push!(bv, nmask[i])
    end
    return bv
end

Base.length(nmask::NMask) = nmask.len

divrem64(i) = i >> 6, i & 0b111111

function block_bit(i)
    d, r = divrem64(Int(i - 1))
    return d + 1, r + 1
end

function block_start(blockid)
    return (blockid - 1) << 6 + 1
end

@inline function Base.getindex(nmask::NMask, i::Integer)
    blockid, bitid = block_bit(i)
    #if !nmask.blockmask[blockid]
    if !IndexableBitVectors.unsafe_getindex(nmask.blockmask, blockid)
        return false
    end
    block = nmask.blocks[IndexableBitVectors.rank1(nmask.blockmask, blockid)]
    return ((block >> (bitid - 1)) & 1) == 1
end

function findnextn(nmask::NMask, i::Integer)
    if i > length(nmask)
        return 0
    end
    blockid, bitid = block_bit(i)
    if nmask.blockmask[blockid]
        # try to find in the current block
        block = nmask.blocks[IndexableBitVectors.rank1(nmask.blockmask, blockid)]
        d = findnext_in_block(block, bitid)
        if d > 0
            # found in the block
            return 64 * (blockid - 1) + d
        end
    end
    # search in the following blocks
    blockid = IndexableBitVectors.search1(nmask.blockmask, blockid + 1)
    if blockid == 0
        return 0
    end
    block = nmask.blocks[IndexableBitVectors.rank1(nmask.blockmask, blockid)]
    d = findnext_in_block(block, 1)
    @assert d > 0
    return 64 * (blockid - 1) + d
end

function findnext_in_block(block, bitid)
    block = block >> (bitid - 1)
    return block == 0 ? 0 : bitid + trailing_zeros(block)
end

function findprevn(nmask::NMask, i::Integer)
    if i ≤ 0
        return 0
    end
    blockid, bitid = block_bit(i)
    if nmask.blockmask[blockid]
        # try to find in the current block
        block = nmask.blocks[IndexableBitVectors.rank1(nmask.blockmask, blockid)]
        d = findprev_in_block(block, bitid)
        if d > 0
            # found in the block
            return 64 * (blockid - 1) + d
        end
    end
    # search in the following blocks
    blockid = IndexableBitVectors.rsearch1(nmask.blockmask, blockid - 1)
    if blockid == 0
        return 0
    end
    block = nmask.blocks[IndexableBitVectors.rank1(nmask.blockmask, blockid)]
    d = findprev_in_block(block, 64)
    @assert d > 0
    return 64 * (blockid - 1) + d
end

function findprev_in_block(block, bitid)
    block = block << (64 - bitid)
    return block == 0 ? 0 : bitid - leading_zeros(block)
end

function hasn_within(nmask::NMask, r::UnitRange{Int})
    # fast heuristic to check that there is no 'N' within the range
    lo, _ = block_bit(first(r))
    hi, _ = block_bit(last(r))
    if (lo     == hi && !nmask.blockmask[lo]) ||
       (lo + 1 == hi && !nmask.blockmask[lo] && !nmask.blockmask[hi])
        return false
    end
    # slower but exact algorithm to check
    return findnextn(nmask, first(r)) in r
end

function make_nbitmap(nmask::NMask, r::UnitRange{Int})
    ns = falses(length(r))
    i = findnextn(nmask, first(r))
    while i in r
        ns[i-first(r)+1] = true
        i = findnextn(nmask, i + 1)
    end
    return ns
end
