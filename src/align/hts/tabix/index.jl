# Tabix Index
# ===========

type Tabix
    format::Int32
    columns::NTuple{3,Int}
    meta::Char
    skip::Int
    names::Vector{String}
    indexes::Vector{Tuple{BinIndex,LinearIndex,Nullable{PseudoBin}}}
    n_no_coor::Nullable{UInt64}
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
    read(filename::AbstractString, ::Type{Tabix})
    read(input::IO, ::Type{Tabix})

Load a Tabix index from a file.
"""
function Tabix(filename::AbstractString)
    return read(filename, Tabix)
end

function Tabix(input::IO)
    return read(input, Tabix)
end

# The BED rule is half-closed-half-open and 0-based like Python.
function is_bed_rule(format)
    return format & 0x10000 != 0
end

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
