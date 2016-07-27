# Hamming Code
# ============
#
# Demultiplexer based on Hamming code.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type HammingCode{k} <: AbstractCode{k}
    codewords::Vector{DNAKmer{k}}
    index::Dict{DNAKmer{k},Int}
end

function HammingCode{k}(codewords::AbstractVector{DNAKmer{k}})
    index = Dict{DNAKmer{k},Int}()
    for (i, codeword) in enumerate(codewords)
        if haskey(index, codeword)
            error("duplicated DNA barcode: ", codeword)
        end
        index[codeword] = i
    end
    return HammingCode{k}(codewords, index)
end

function codewords(code::HammingCode)
    return code.codewords
end

function demultiplex{k,A<:DNAAlphabet}(code::HammingCode{k}, seq::BioSequence{A})
    barcode = extractbarcode(DNAKmer{k}, seq, 1)
    i = get(code.index, barcode, 0)
    if i != 0
        return i
    end
    i_min = 0
    dist_min = typemax(Int)
    for (i, codeword) in enumerate(code.codewords)
        dist = mismatches(codeword, barcode)
        if dist < dist_min
            i_min = i
            dist_min = dist
        end
    end
    return i_min
end

function extractbarcode{k}(::Type{DNAKmer{k}}, seq, from)
    x::UInt64 = 0
    for i in from:from+k-1
        x <<= 2
        nt = seq[i]
        if isambiguous(nt)
            # fill ambiguous nucleotides with 'A's
            x |= UInt64(DNA_A)
        else
            x |= UInt64(nt)
        end
    end
    return DNAKmer{k}(x)
end
