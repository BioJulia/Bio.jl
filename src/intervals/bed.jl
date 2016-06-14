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

@compat function Base.:(==)(a::BEDMetadata, b::BEDMetadata)
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

function Base.open(filepath::AbstractString, mode::AbstractString, ::Type{BED};
                   n_fields::Integer=-1)
    io = open(filepath, mode)
    if mode[1] == 'r'
        return open(BufferedInputStream(io), BED)
    elseif mode[1] ∈ ('w', 'o')
        return BEDWriter(io, n_fields)
    end
    error("invalid open mode")
end


# Parser
# ------

type BEDParser <: AbstractParser
    state::Ragel.State

    # intermediate values used during parsing
    red::Float32
    green::Float32
    blue::Float32
    block_size_idx::Int
    block_first_idx::Int

    function BEDParser(input::BufferedInputStream)
        return new(Ragel.State(bedparser_start, input), 0, 0, 0, 1, 1)
    end
end

function Intervals.metadatatype(::BEDParser)
    return BEDMetadata
end

function Base.eltype(::Type{BEDParser})
    return BEDInterval
end

function Base.eof(parser::BEDParser)
    return eof(parser.state.stream)
end

function Base.open(input::BufferedInputStream, ::Type{BED})
    return BEDParser(input)
end

function IntervalCollection(interval_stream::BEDParser)
    intervals = collect(BEDInterval, interval_stream)
    return IntervalCollection{BEDMetadata}(intervals, true)
end

include("bed-parser.jl")

# Writer
# ------

type BEDWriter{T<:IO} <: AbstractWriter
    output::T

    # number of written additional fields (-1: inferred from the first interval)
    n_fields::Int
end

function Base.flush(writer::BEDWriter)
    # TODO: This can be removed on Julia v0.5
    # (because flush will be defined for IOBuffer).
    if applicable(flush, writer.output)
        flush(writer.output)
    end
end

Base.close(writer::BEDWriter) = close(writer.output)

function Base.write(writer::BEDWriter, interval::BEDInterval)
    if writer.n_fields == -1
        writer.n_fields = interval.metadata.used_fields
    end

    if interval.metadata.used_fields < writer.n_fields
        error("interval doesn't have enough additional fields")
    end

    output = writer.output
    n_fields = writer.n_fields
    n = 0
    n += write(
        output,
        interval.seqname, '\t',
        dec(interval.first - 1), '\t',
        dec(interval.last))

    if n_fields ≥ 1
        n += write(output, '\t', interval.metadata.name)
    end

    if n_fields ≥ 2
        n += write(output, '\t', dec(interval.metadata.score))
    end

    if n_fields ≥ 3
        n += write(output, '\t', Char(interval.strand))
    end

    if n_fields ≥ 4
        n += write(output, '\t', dec(interval.metadata.thick_first - 1))
    end

    if n_fields ≥ 5
        n += write(output, '\t', dec(interval.metadata.thick_last))
    end

    if n_fields ≥ 6
        rgb = interval.metadata.item_rgb
        n += write(
            output, '\t',
            dec(round(Int, 255 * rgb.r)), ',',
            dec(round(Int, 255 * rgb.g)), ',',
            dec(round(Int, 255 * rgb.b)))
    end

    if n_fields ≥ 7
        n += write(output, '\t', dec(interval.metadata.block_count))
    end

    if n_fields ≥ 8
        block_sizes = interval.metadata.block_sizes
        if !isempty(block_sizes)
            n += write(output, '\t', dec(block_sizes[1]))
            for i in 2:length(block_sizes)
                n += write(output, ',', dec(block_sizes[i]))
            end
        end
    end

    if n_fields ≥ 9
        block_firsts = interval.metadata.block_firsts
        if !isempty(block_firsts)
            n += write(output, '\t', dec(block_firsts[1] - 1))
            for i in 2:length(block_firsts)
                n += write(output, ',', dec(block_firsts[i] - 1))
            end
        end
    end

    n += write(output, '\n')
    return n
end

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
