# CIGAR
# =====
#
# CIGAR string.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type CIGAR
    data::Vector{UInt64}
end

# empty CIGAR string
function CIGAR()
    return CIGAR(UInt64[])
end

# run-length encoding of operations
function CIGAR(oplens::Vector{Tuple{Operation,Integer}})
    data = Vector{UInt64}(length(oplens))
    for (op, len) in oplens
        data[i] = encode_cigar(op, len)
    end
    return CIGAR(data)
end

# parse CIGAR string
function CIGAR(s::AbstractString)
    data = UInt64[]
    len = 0
    for c in s
        if isnumber(c)
            len = 10len + (c - '0')
        else
            op = convert(Operation, c)
            push!(data, encode_cigar(op, len))
            len = 0
        end
    end
    return CIGAR(data)
end

function Base.getindex(cigar::CIGAR, i::Integer)
    x = cigar.data[i]
    return cigarop(x), cigarlen(x)
end

function Base.push!(cigar::CIGAR, oplen::Tuple{Operation,Integer})
    push!(cigar.data, encode_cigar(oplen[1], oplen[2]))
    return cigar
end

Base.eltype(::Type{CIGAR}) = Tuple{Operation,Int}
Base.length(cigar::CIGAR) = length(cigar.data)
Base.start(::CIGAR) = 1
Base.done(cigar::CIGAR, i) = i > length(cigar)
Base.next(cigar::CIGAR, i) = cigar[i], i + 1

function Base.show(io::IO, cigar::CIGAR)
    print(io, "CIGAR(\"", string(cigar), "\")")
end

function Base.print(io::IO, cigar::CIGAR)
    for (op, len) in cigar
        print(io, len, op)
    end
end

function encode_cigar(op::Operation, len::Integer)
    if !isvalid(op)
        throw(ArgumentError("invalid operation"))
    elseif len > 2^60
        throw(ArgumentError("too long operation length"))
    end
    return UInt64(len) << 4 | UInt64(op)
end
cigarlen(x::UInt64) = Int(x >> 4)
cigarop(x::UInt64) = Operation(x & 0x0f)
