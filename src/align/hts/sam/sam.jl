# SAM
# ===

"""
The SAM file format.
"""
immutable SAM <: Bio.IO.FileFormat end

function checkkeytag(tag)
    if !isa(tag, AbstractString)
        throw(ArgumentError("tag must be a string"))
    elseif length(tag) != 2
        throw(ArgumentError("tag length must be 2"))
    end
    return
end

const auxtype = Dict{UInt8,DataType}(
    UInt8('A') => Char,
    UInt8('c') => Int8,
    UInt8('C') => UInt8,
    UInt8('s') => Int16,
    UInt8('S') => UInt16,
    UInt8('i') => Int32,
    UInt8('I') => UInt32,
    UInt8('f') => Float32,
    UInt8('d') => Float64,
    UInt8('Z') => String)

const auxtypechar = Dict{DataType,Char}(
    Char    => 'A',
    Int8    => 'c',
    UInt8   => 'C',
    Int16   => 's',
    UInt16  => 'S',
    Int32   => 'i',
    UInt32  => 'I',
    Float32 => 'f',
    Float64 => 'd',
    String  => 'Z')

# Parse hexadecimal byte format string.
function parse_hexbytearray(s)
    return [parse(UInt8, s[i:i+1], 16) for i in 1:2:endof(s)]
end

function parse_keyvals(s)
    ret = Dict{String,Any}()
    for pair in split(s, '\t')
        if pair[3] != ':'
            error("':' is expected")
        end
        ret[pair[1:2]] = pair[4:end]
    end
    return ret
end

include("flags.jl")
include("header.jl")
include("record.jl")
include("reader.jl")
include("parser.jl")
include("writer.jl")
