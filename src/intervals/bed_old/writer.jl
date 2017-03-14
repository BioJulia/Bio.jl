# BED Writer
# ==========

"""
    BEDWriter(output::IO)

Create a data writer of the BED file format.

# Arguments
* `output`: data sink
"""
type BEDWriter{T<:IO} <: Bio.IO.AbstractWriter
    output::T

    # number of written additional fields (-1: inferred from the first interval)
    n_fields::Int
end

function Bio.IO.stream(writer::BEDWriter)
    return writer.output
end

function BEDWriter(output::IO, n_fields::Integer=-1)
    return BEDWriter{typeof(output)}(output, n_fields)
end

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
