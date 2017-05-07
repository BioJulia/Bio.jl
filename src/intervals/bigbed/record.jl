# BigBed Record
# =============

# type Record is defined in reader.jl.

"""
    BigBed.Record()

Create an unfilled bigBed record.
"""
function Record()
    return Record(
        # chromid, chromstart, chromend
        0, 0, 0,
        # data, filled range, #columns
        UInt8[], 1:0, 0,
        # indexes
        1:0, 1:0, 0,
        1:0, 1:0, 1:0,
        1:0, UnitRange{Int}[], UnitRange{Int}[])
end

function Base.convert(::Type{Bio.Intervals.Interval}, record::Record)
    s = hasstrand(record) ? strand(record) : Bio.Intervals.STRAND_BOTH
    return Bio.Intervals.Interval(chrom(record), chromstart(record), chromend(record), s, record)
end

function Base.convert(::Type{Bio.Intervals.Interval{Record}}, record::Record)
    return convert(Bio.Intervals.Interval, record)
end

function Base.copy(record::Record)
    copy = Record(
        record.chromid,
        record.chromstart,
        record.chromend,
        record.data[record.filled],
        record.filled,
        record.ncols,
        record.name,
        record.score,
        record.strand,
        record.thickstart,
        record.thickend,
        record.itemrgb,
        record.blockcount,
        Base.copy(record.blocksizes),
        Base.copy(record.blockstarts))
    if isdefined(record, :reader)
        copy.reader = record.reader
    end
    return copy
end

function initialize!(record::Record)
    record.chromid = 0
    record.chromstart = 0
    record.chromend = 0
    record.filled = 1:0
    record.ncols = 0
    record.name = 1:0
    record.score = 1:0
    record.strand = 0
    record.thickstart = 1:0
    record.thickend = 1:0
    record.itemrgb = 1:0
    record.blockcount = 1:0
    empty!(record.blocksizes)
    empty!(record.blockstarts)
    return record
end

function Bio.isfilled(record::Record)
    return !isempty(record.filled)
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   chromosome: ", chrom(record))
        println(io, "        start: ", chromstart(record))
          print(io, "          end: ", chromend(record))
        if hasname(record)
            println(io)
            print(io, "         name: ", name(record))
        end
        if hasscore(record)
            println(io)
            print(io, "        score: ", score(record))
        end
        if hasstrand(record)
            println(io)
            print(io, "       strand: ", strand(record))
        end
        if hasthickstart(record)
            println(io)
            print(io, "  thick start: ", thickstart(record))
        end
        if hasthickend(record)
            println(io)
            print(io, "    thick end: ", thickend(record))
        end
        if hasitemrgb(record)
            println(io)
            print(io, "     item RGB: ", itemrgb(record))
        end
        if hasblockcount(record)
            bcount = blockcount(record)
            println(io)
            print(io, "  block count: ", bcount)
            for (i, (bsize, bstart)) in enumerate(zip(blocksizes(record), blockstarts(record)))
                println(io)
                print(io, "      [$i]: size=$(bsize), start=$(bstart)")
            end
        end
    else
        print(io, " <not filled>")
    end
end


# Accessor functions
# ------------------

"""
    chromid(record::Record)::UInt32

Get the chromosome ID of `record`.
"""
function chromid(record::Record)::UInt32
    checkfilled(record)
    return record.chromid
end

"""
    chrom(record::Record)::String

Get the chromosome name of `record`.
"""
function chrom(record::Record)::String
    checkfilled(record)
    return record.reader.chrom_names[chromid(record)]
end

"""
    chromstart(record::Record)::Int

Get the start position of `record`.
"""
function chromstart(record::Record)::Int
    checkfilled(record)
    return record.chromstart + 1
end

"""
    chromend(record::Record)::Int

Get the end position of `record`.
"""
function chromend(record::Record)::Int
    checkfilled(record)
    return record.chromend % Int
end

"""
    name(record::Record)::String

Get the name of `record`.
"""
function name(record::Record)::String
    checkfilled(record)
    if !hasname(record)
        missingerror(:name)
    end
    return String(record.data[record.name])
end

function hasname(record::Record)
    return record.ncols ≥ 4
end

"""
    score(record::Record)::Int

Get the score between 0 and 1000.
"""
function score(record::Record)::Int
    checkfilled(record)
    if !hasscore(record)
        missingerror(:score)
    end
    return BED.unsafe_parse_decimal(Int, record.data, record.score)
end

function hasscore(record::Record)
    return record.ncols ≥ 5
end

"""
    strand(record::Record)::Bio.Intervals.Strand

Get the strand of `record`.
"""
function strand(record::Record)::Bio.Intervals.Strand
    checkfilled(record)
    if !hasstrand(record)
        missingerror(:strand)
    end
    return convert(Bio.Intervals.Strand, Char(record.data[record.strand]))
end

function hasstrand(record::Record)
    return record.ncols ≥ 6
end

function Bio.Intervals.strand(record::Record)
    return strand(record)
end

"""
    thickstart(record::Record)::Int

Get the starting position at which `record` is drawn thickly.

Note that the first base is numbered 1.
"""
function thickstart(record::Record)::Int
    checkfilled(record)
    if !hasthickstart(record)
        missingerror(:thickstart)
    end
    return BED.unsafe_parse_decimal(Int, record.data, record.thickstart) + 1
end

function hasthickstart(record::Record)
    return record.ncols ≥ 7
end

"""
    thickend(record::Record)::Int

Get the end position at which `record` is drawn thickly.
"""
function thickend(record::Record)::Int
    checkfilled(record)
    if !hasthickend(record)
        missingerror(:thickend)
    end
    return BED.unsafe_parse_decimal(Int, record.data, record.thickend)
end

function hasthickend(record::Record)
    return record.ncols ≥ 8
end

"""
    itemrgb(record::Record)::ColorTypes.RGB

Get the RGB value of `record`.

The return type is defined in [ColorTypes.jl](https://github.com/JuliaGraphics/ColorTypes.jl).
"""
function itemrgb(record::Record)::ColorTypes.RGB
    checkfilled(record)
    if !hasitemrgb(record)
        missingerror(:itemrgb)
    end
    return BED.parse_rgbcolor(record.data, record.itemrgb)
end

function hasitemrgb(record::Record)
    return record.ncols ≥ 9
end

"""
    blockcount(record::Record)::Int

Get the number of blocks (exons) in `record`.
"""
function blockcount(record::Record)::Int
    checkfilled(record)
    if !hasblockcount(record)
        missingerror(:blockcount)
    end
    return BED.unsafe_parse_decimal(Int, record.data, record.blockcount)
end

function hasblockcount(record::Record)
    return record.ncols ≥ 10
end

"""
    blocksizes(record::Record)::Vector{Int}

Get the block (exon) sizes of `record`.
"""
function blocksizes(record::Record)::Vector{Int}
    checkfilled(record)
    if !hasblocksizes(record)
        missingerror(:blocksizes)
    end
    return [BED.unsafe_parse_decimal(Int, record.data, r) for r in record.blocksizes]
end

function hasblocksizes(record::Record)
    return record.ncols ≥ 11
end

"""
    blockstarts(record::Record)::Vector{Int}

Get the block (exon) starts of `record`.

Note that the first base is numbered 1.
"""
function blockstarts(record::Record)::Vector{Int}
    checkfilled(record)
    if !hasblockstarts(record)
        missingerror(:blockstarts)
    end
    return [BED.unsafe_parse_decimal(Int, record.data, r) + 1 for r in record.blockstarts]
end

function hasblockstarts(record::Record)
    return record.ncols ≥ 12
end

"""
    optionals(record::Record)::Vector{String}

Get optional fields as strings.
"""
function optionals(record::Record)::Vector{String}
    checkfilled(record)
    ret = String[]
    p = 13
    p_end = search(record.data, '\0', p)
    while p < p_end
        p_delim = find(record.data, '\t', p)
        push!(ret, String(record.data[p:p_delim-1]))
        p = p_delim + 1
    end
    return ret
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled bigBed record"))
    end
end
