# Demultiplex
# ===========
#
# Demultiplexer for multiplexed sequencing.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract AbstractCode{k} <: AbstractVector{DNAKmer{k}}

"""
    codewords(code::AbstractCode)

Return codewords (or DNA barcodes) used for demultiplexing sequences.
"""
function codewords end

"""
    demultiplex(code::AbstractCode, seq::BioSequence)

Demultiplex `seq` based on the prefix `seq` and return the barcode index if it
succeeds; otherwise return `0`.
"""
function demultiplex end

function Base.getindex(code::AbstractCode, i::Integer)
    return codewords(code)[i]
end

function Base.size(code::AbstractCode)
    return size(codewords(code))
end

include("hamming.jl")
include("sequence_levenshtein.jl")
