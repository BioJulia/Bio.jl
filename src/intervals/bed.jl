# BED
# ===
#
# Reader and writer of the BED file format.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable BED <: FileFormat end

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

function Base.(:(==))(a::BEDMetadata, b::BEDMetadata)
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

# TODO
#function show(io::IO, metadata::BEDMetadata)
#end

"An `Interval` with associated metadata from a BED file"
typealias BEDInterval Interval{BEDMetadata}

"""
Write a BEDInterval in BED format.
"""
function Base.write(out::IO, interval::BEDInterval)
    print(out, interval.seqname, '\t', interval.first - 1, '\t', interval.last)
    write_optional_fields(out, interval)
    write(out, '\n')
end

function write_optional_fields(out::IO, interval::BEDInterval, leadingtab::Bool=true)
    if interval.metadata.used_fields >= 1
        if leadingtab
            write(out, '\t')
        end
        print(out, interval.metadata.name)
    else return end

    if interval.metadata.used_fields >= 2
        print(out, '\t', interval.metadata.score)
    else return end

    if interval.metadata.used_fields >= 3
        print(out, '\t', interval.strand)
    else return end

    if interval.metadata.used_fields >= 4
        print(out, '\t', interval.metadata.thick_first - 1)
    else return end

    if interval.metadata.used_fields >= 5
        print(out, '\t', interval.metadata.thick_last)
    else return end

    if interval.metadata.used_fields >= 6
        item_rgb = interval.metadata.item_rgb
        print(out, '\t',
              round(Int, 255 * item_rgb.r), ',',
              round(Int, 255 * item_rgb.g), ',',
              round(Int, 255 * item_rgb.b))
    else return end

    if interval.metadata.used_fields >= 7
        print(out, '\t', interval.metadata.block_count)
    else return end

    if interval.metadata.used_fields >= 8
        block_sizes = interval.metadata.block_sizes
        if !isempty(block_sizes)
            print(out, '\t', block_sizes[1])
            for i in 2:length(block_sizes)
                print(out, ',', block_sizes[i])
            end
        end
    else return end

    if interval.metadata.used_fields >= 9
        block_firsts = interval.metadata.block_firsts
        if !isempty(block_firsts)
            print(out, '\t', block_firsts[1] - 1)
            for i in 2:length(block_firsts)
                print(out, ',', block_firsts[i] - 1)
            end
        end
    end
end
