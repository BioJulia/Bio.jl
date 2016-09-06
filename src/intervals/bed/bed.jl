# BED
# ===
#
# Reader and writer of the BED file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""Metadata for BED interval records"""
type BEDMetadata
    used_fields::Int # how many of the first n fields are used
    name::StringField
    score::Int
    thick_first::Int
    thick_last::Int
    item_rgb::RGB{Float32}
    block_count::Int
    block_sizes::Vector{Int}
    block_firsts::Vector{Int}
end

function BEDMetadata()
    return BEDMetadata(0, StringField(), 0, 0, 0, RGB{Float32}(0.0, 0.0, 0.0),
                       0, Int[], Int[])
end

function Base.copy(metadata::BEDMetadata)
    return BEDMetadata(
        metadata.used_fields, copy(metadata.name),
        metadata.score, metadata.thick_first, metadata.thick_last,
        metadata.item_rgb, metadata.block_count,
        metadata.block_sizes[1:metadata.block_count],
        metadata.block_firsts[1:metadata.block_count])
end

function Base.:(==)(a::BEDMetadata, b::BEDMetadata)
    if a.used_fields != b.used_fields
        return false
    end

    n = a.used_fields
    ans = (n < 1 || a.name == b.name) &&
          (n < 2 || a.score == b.score) &&
          (n < 4 || a.thick_first == b.thick_first) &&
          (n < 5 || a.thick_last == b.thick_last) &&
          (n < 6 || a.item_rgb == b.item_rgb) &&
          (n < 7 || a.block_count == b.block_count)
    if !ans
        return false
    end

    if n >= 8
        for i in 1:a.block_count
            if a.block_sizes[i] != b.block_sizes[i]
                return false
            end
        end
    end

    if n >= 9
        for i in 1:a.block_count
            if a.block_sizes[i] != b.block_sizes[i]
                return false
            end
        end
    end

    return true
end

"An `Interval` with associated metadata from a BED file"
typealias BEDInterval Interval{BEDMetadata}

include("reader.jl")
include("parser.jl")
include("writer.jl")
