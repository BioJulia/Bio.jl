# BAI
# ===
#
# Index for BAM files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# An index type for the BAM file format.
type BAI
    # BGZF file index
    index::Bio.Intervals.BGZFIndex

    # number of unmapped reads
    n_no_coors::Nullable{Int}
end

"""
    BAI(filename::AbstractString)

Load a BAI index from `filename`.
"""
function BAI(filename::AbstractString)
    return open(read_bai, filename)
end

"""
    BAI(input::IO)

Load a BAI index from `input`.
"""
function BAI(input::IO)
    return read_bai(input)
end

# Read a BAI object from `input`.
function read_bai(input::IO)
    # check magic bytes
    B = read(input, UInt8)
    A = read(input, UInt8)
    I = read(input, UInt8)
    x = read(input, UInt8)
    if B != UInt8('B') || A != UInt8('A') || I != UInt8('I') || x != 0x01
        error("input is not a valid BAI file")
    end

    # read contents
    n_refs = read(input, Int32)
    index = Bio.Intervals.read_bgzfindex(input, n_refs)
    if !eof(input)
        n_no_coors = Nullable{Int}(read(input, UInt64))
    else
        n_no_coors = Nullable{Int}()
    end

    return BAI(index, n_no_coors)
end
