# FASTA Index
# ===========
#
# Index for random access to FASTA files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# http://www.htslib.org/doc/faidx.html
type FASTAIndex
    names::Vector{String}
    lengths::Vector{Int}
    offsets::Vector{Int}
    linebases::Vector{Int}
    linewidths::Vector{Int}
end

function FASTAIndex(filepath::AbstractString)
    return open(read_faidx, filepath)
end

function read_faidx(input::IO)
    names = String[]
    lengths = Int[]
    offsets = Int[]
    linebases = Int[]
    linewidths = Int[]
    for line in eachline(input)
        values = split(chomp(line), '\t')
        name = values[1]
        length = parse(Int, values[2])
        offset = parse(Int, values[3])
        linebase = parse(Int, values[4])
        linewidth = parse(Int, values[5])
        push!(names, name)
        push!(lengths, length)
        push!(offsets, offset)
        push!(linebases, linebase)
        push!(linewidths, linewidth)
    end
    return FASTAIndex(names, lengths, offsets, linebases, linewidths)
end

# Set the reading position of `input` to the starting position of the record `name`.
function seekrecord(input::IO, fai::FASTAIndex, name::AbstractString)
    i = findfirst(fai.names, name)
    if i == 0
        error("sequence \"", name, "\" is not in the index")
    end
    offset = fai.offsets[i]
    n_back = 100
    @label seekback
    seek(input, max(offset - n_back, 0))
    data = UInt8[]
    while position(input) < offset
        push!(data, read(input, UInt8))
    end
    for j in endof(data):-1:1
        if data[j] == UInt8('>') && (offset ≤ n_back || (j ≥ 2 && data[j-1] == UInt8('\n')))
            seek(input, offset - (endof(data) - j + 1))
            return
        end
    end
    if n_back ≥ offset
        # reached the starting position of the input
        error("failed to find the starting position of the record")
    end
    n_back *= 2
    @goto seekback
end

