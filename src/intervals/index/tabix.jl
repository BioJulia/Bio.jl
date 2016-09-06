# Tabix
# =====
#
# Generic index for tab-delimited files.
#
# Li, Heng. "Tabix: fast retrieval of sequence features from generic
# TAB-delimited files." Bioinformatics 27.5 (2011): 718-719.
# Specification: http://samtools.github.io/hts-specs/tabix.pdf
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# An index type for tab-delimited files.
type Tabix
    # file format
    #   * 0: generic
    #   * 1: SAM
    #   * 2: VCF
    # note: `format & 0x10000 != 0` indicates the BED rule.
    format::Int32

    # triplet of columns (sequence name, start of a region, end of a region)
    columns::NTuple{3,Int}

    # leading character for comment lines
    meta::Char

    # number of lines to skip at the beginning
    skip::Int

    # sequence names
    names::Vector{String}

    # BGZF file index
    index::BGZFIndex

    # number of unmapped reads
    n_no_coor::Nullable{Int}
end

function Base.show(io::IO, index::Tabix)
    println(io, summary(index), ":")
    println(io, "  format: ", format2str(index.format))
    println(io, "  columns: ", index.columns)
    println(io, "  meta char: '", index.meta, "'")
    println(io, "  skip lines: ", index.skip)
      print(io, "  names: ", index.names)
end

"""
    Tabix(filename::AbstractString)
    Tabix(input::IO)

Load a Tabix index from `filename` or `input`.
"""
function Tabix(filename::AbstractString)
    return open(read_tabix, filename)
end

function Tabix(input::IO)
    return read_tabix(input)
end

"""
    overlapchunks(tabix::Tabix, interval::Interval)

Return chunks possibly overlapping with the range specified by `interval`.

Note that records within the returned chunks are not guaranteed to actually
overlap the query interval.
"""
function overlapchunks(tabix::Tabix, interval::Interval)
    seqid = findfirst(tabix.names, seqname(interval))
    if seqid == 0
        throw(ArgumentError("sequence name $(seqname) is not included in the index"))
    end
    return overlapchunks(tabix.index, seqid, interval.first, interval.last)
end

# Check if `format` follows the BED rule (half-closed-half-open and 0-based).
function is_bed_rule(format)
    return format & 0x10000 != 0
end

# Convert a tabix file format integer to a string.
function format2str(format)
    if format == 1
        return "SAM"
    elseif format == 2
        return "VCF"
    else
        if is_bed_rule(format)
            return "generic (BED rule)"
        else
            return "generic"
        end
    end
end

# Read a Tabix object from `input_`.
function read_tabix(input_::IO)
    input = BGZFStream(input_)

    # check magic bytes
    T = read(input, UInt8)
    B = read(input, UInt8)
    I = read(input, UInt8)
    x = read(input, UInt8)
    if T != UInt8('T') || B != UInt8('B') || I != UInt8('I') || x != 0x01
        error("invalid tabix magic bytes")
    end

    # read contents
    n_refs = read(input, Int32)
    format = read(input, Int32)
    col_seq = read(input, Int32)
    col_beg = read(input, Int32)
    col_end = read(input, Int32)
    meta = read(input, Int32)
    skip = read(input, Int32)
    l_nm = read(input, Int32)
    data = read(input, UInt8, l_nm)
    names = split(String(data), '\0', keep=false)
    if length(names) != n_refs
        error("the number of sequence names doesn't match the expacted value")
    end
    index = read_bgzfindex(input, n_refs)
    if !eof(input)
        n_no_coor = Nullable{Int}(read(input, UInt64))
    else
        n_no_coor = Nullable{Int}()
    end

    return Tabix(
        format,
        (col_seq, col_beg, col_end),
        meta,
        skip,
        names,
        index,
        n_no_coor)
end
