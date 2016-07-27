# Sequence-Levenshtein Code
# =========================
#
# Demultiplexer based on Sequence-Levenshtein code.
#
# Buschmann, Tilo, and Leonid V. Bystrykh. "Levenshtein error-correcting
# barcodes for multiplexed DNA sequencing." BMC bioinformatics 14.1 (2013): 272.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type SequenceLevenshteinCode{k} <: AbstractCode{k}
    codewords::Vector{DNAKmer{k}}
    index::Dict{DNAKmer{k},Int}
end

function SequenceLevenshteinCode{k}(codewords::AbstractVector{DNAKmer{k}})
    index = Dict{DNAKmer{k},Int}()
    for (i, codeword) in enumerate(codewords)
        if haskey(index, codeword)
            error("duplicated DNA barcode: ", codeword)
        end
        index[codeword] = i
    end
    return SequenceLevenshteinCode{k}(codewords, index)
end

function codewords(code::SequenceLevenshteinCode)
    return code.codewords
end

function demultiplex{k,A<:DNAAlphabet}(code::SequenceLevenshteinCode{k}, seq::BioSequence{A})
    barcode = extractbarcode(DNAKmer{k}, seq, 1)
    i = get(code.index, barcode, 0)
    if i != 0
        # perfect match
        return i
    end
    # distance matrix used for dynamic programming
    distmtx = Matrix{Int}(k + 1, length(seq) + 1)
    i_min = 0
    dist_min = typemax(Int)
    for (i, codeword) in enumerate(code.codewords)
        dist = seqlevdist!(distmtx, codeword, seq)
        if dist < dist_min
            i_min = i
            dist_min = dist
        end
    end
    return i_min
end

# Compute Sequence-Levenshtein distance.
function seqlevdist!(distmtx, codeword, seq)
    k = length(codeword)
    l = length(seq)
    @assert size(distmtx) == (k + 1, l + 1)

    for i in 0:k
        distmtx[i+1,1] = i
    end
    for j in 1:l
        distmtx[1,j+1] = j
        for i in 1:k
            distmtx[i+1,j+1] = min(
                distmtx[i,j+1] + 1,
                distmtx[i+1,j] + 1,
                distmtx[i,j]   + ifelse(codeword[i] == seq[j], 0, 1))
        end
    end

    mindist = distmtx[end,end]
    for i in 0:k
        mindist = min(mindist, distmtx[i+1,end])
    end
    for j in 0:l
        mindist = min(mindist, distmtx[end,j+1])
    end

    return mindist
end
