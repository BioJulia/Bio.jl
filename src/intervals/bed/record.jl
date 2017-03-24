# BED Record
# ==========

type Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # number of columns
    ncols::Int
    # indexes
    chrom::UnitRange{Int}
    chromstart::UnitRange{Int}
    chromend::UnitRange{Int}
    name::UnitRange{Int}
    score::UnitRange{Int}
    strand::Int
    thickstart::UnitRange{Int}
    thickend::UnitRange{Int}
    itemrgb::UnitRange{Int}
    blockcount::UnitRange{Int}
    blocksizes::Vector{UnitRange{Int}}
    blockstarts::Vector{UnitRange{Int}}
end

"""
    BED.Record()

Create an unfilled BED record.
"""
function Record()
    return Record(
        UInt8[], 1:0, 0,
        # chrom-score
        1:0, 1:0, 1:0, 1:0, 1:0,
        # strand-itemrgb
        0, 1:0, 1:0, 1:0,
        # blockcount-blockstarts
        1:0, UnitRange{Int}[], UnitRange{Int}[])
end

"""
    BED.Record(data::Vector{UInt8})

Create a BED record object from `data`.

This function verifies and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(
        data, 1:0, 0,
        # chrom-score
        1:0, 1:0, 1:0, 1:0, 1:0,
        # strand-itemrgb
        0, 1:0, 1:0, 1:0,
        # blockcount-blockstarts
        1:0, UnitRange{Int}[], UnitRange{Int}[])
    index!(record)
    return record
end

"""
    BED.Record(str::AbstractString)

Create a BED record object from `str`.

This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return convert(Record, convert(Vector{UInt8}, str))
end

function initialize!(record::Record)
    record.filled = 1:0
    record.ncols = 0
    record.chrom = 1:0
    record.chromstart = 1:0
    record.chromend = 1:0
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

function Base.convert(::Type{Bio.Intervals.Interval}, record::Record)
    name = Bio.seqname(record)
    lpos = Bio.leftposition(record)
    rpos = Bio.rightposition(record)
    strd = hasstrand(record) ? Bio.Intervals.strand(record) : Bio.Intervals.STRAND_BOTH
    return Bio.Intervals.Interval(name, lpos, rpos, strd, record)
end

function Base.convert(::Type{Bio.Intervals.Interval{Record}}, record::Record)
    return convert(Bio.Intervals.Interval, record)
end

function isfilled(record::Record)
    return !isempty(record.filled)
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = record1.filled
        r2 = record2.filled
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.ncols,
        record.chrom,
        record.chromstart,
        record.chromend,
        record.name,
        record.score,
        record.strand,
        record.thickstart,
        record.thickend,
        record.itemrgb,
        record.blockcount,
        copy(record.blocksizes),
        copy(record.blockstarts))
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
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
    chrom(record::Record)::String

Get the chromosome name of `record`.
"""
function chrom(record::Record)::String
    checkfilled(record)
    return String(record.data[record.chrom])
end

function haschrom(record::Record)
    return isfilled(record)
end

function Bio.seqname(record::Record)
    return chrom(record)
end

function Bio.hasseqname(record::Record)
    return haschrom(record)
end

"""
    chromstart(record::Record)::Int

Get the starting position of `record`.

Note that the first base is numbered 1.
"""
function chromstart(record::Record)::Int
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.chromstart) + 1
end

function haschromstart(record::Record)
    return isfilled(record)
end

function Bio.leftposition(record::Record)
    return chromstart(record)
end

function Bio.hasleftposition(record::Record)
    return haschromstart(record)
end

"""
    chromend(record::Record)::Int

Get the end position of `record`.
"""
function chromend(record::Record)::Int
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.chromend)
end

function haschromend(record::Record)
    return isfilled(record)
end

function Bio.rightposition(record::Record)
    return chromend(record)
end

function Bio.hasrightposition(record::Record)
    return haschromend(record)
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
    return unsafe_parse_decimal(Int, record.data, record.score)
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
    return unsafe_parse_decimal(Int, record.data, record.thickstart) + 1
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
    return unsafe_parse_decimal(Int, record.data, record.thickend)
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
    return parse_rgbcolor(record.data, record.itemrgb)
end

function hasitemrgb(record::Record)
    return record.ncols ≥ 9
end

function parse_rgbcolor(data::Vector{UInt8}, range::UnitRange{Int})
    function searchcomma(s)
        i = s
        while i ≤ last(range)
            if data[i] == UInt8(',')
                break
            end
            i += 1
        end
        return i
    end
    lo = first(range)
    hi = searchcomma(lo)
    r = unsafe_parse_byte(data, lo:hi-1)
    if hi > last(range)
        # single value
        g = b = r
    else
        # triplet
        lo = hi + 1
        hi = searchcomma(lo)
        g = unsafe_parse_byte(data, lo:hi-1)
        lo = hi + 1
        hi = searchcomma(lo)
        b = unsafe_parse_byte(data, lo:hi-1)
    end
    return ColorTypes.RGB(reinterpret(N0f8, r), reinterpret(N0f8, g), reinterpret(N0f8, b))
end

function unsafe_parse_byte(data::Vector{UInt8}, range::UnitRange{Int})
    val::UInt8 = 0x00
    for i in range
        val = val * 0x0a + (data[i] - UInt8('0'))
    end
    return val
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
    return unsafe_parse_decimal(Int, record.data, record.blockcount)
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
    return [unsafe_parse_decimal(Int, record.data, r) for r in record.blocksizes]
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
    return [unsafe_parse_decimal(Int, record.data, r) + 1 for r in record.blockstarts]
end

function hasblockstarts(record::Record)
    return record.ncols ≥ 12
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled BED record"))
    end
end

# r"[-+]?[0-9]+" must match `data[range]`.
function unsafe_parse_decimal{T<:Signed}(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int})
    lo = first(range)
    if data[lo] == UInt8('-')
        sign = T(-1)
        lo += 1
    elseif data[lo] == UInt8('+')
        sign = T(+1)
        lo += 1
    else
        sign = T(+1)
    end
    x = zero(T)
    @inbounds for i in lo:last(range)
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return sign * x
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end
