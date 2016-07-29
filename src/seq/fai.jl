# FAIndex
# =======
#
# Index for random access to FASTA files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# http://www.htslib.org/doc/faidx.html
type FAIndex
    names::Vector{ASCIIString}
    lengths::Vector{Int}
    offsets::Vector{Int}
    linebases::Vector{Int}
    linewidths::Vector{Int}
end

function FAIndex(filepath::AbstractString)
    return open(read_faidx, filepath)
end

function read_faidx(input::IO)
    names = ASCIIString[]
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
    return FAIndex(names, lengths, offsets, linebases, linewidths)
end

# Set the reading position of `input` to the beginning of a record.
function seekseq(input::IO, fai::FAIndex, name::AbstractString)
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
    reverse!(data)
    if offset == minimum(fai.offsets)
        s = first(search(bytestring(data), r">"))
    else
        s = first(search(bytestring(data), r">\n"))
    end
    if s == 0
        n_back *= 2
        @goto seekback
    end
    seek(input, offset - s)
    return
end

