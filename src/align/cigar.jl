
immutable CIGAR
    OP::Operation
    Size::Int
end

function CIGAR(op::Char, size::Int)
    return CIGAR(Operation(op), size)
end

function convert(::Type{CIGAR}, str::String)
    return CIGAR(str[end], parse(Int, str[1:end-1]))
end

function convert(::Type{String}, cigar::CIGAR)
    return string(cigar.Size, Char(cigar.OP))
end

function show(io::IO, cigar::CIGAR)
    write(io, convert(String, cigar))
end

typealias CIGARString Vector{CIGAR}

function convert(::Type{CIGARString}, str::String)
    matches = matchall(r"(\d+)(\w)", str)
    cigarString = Vector{CIGAR}(length(matches))
    @inbounds for i in 1:length(matches)
        cigarString[i] = CIGAR(matches[i])
    end
    return cigarString
end

macro cigar_str(str)
    return CIGARString(str)
end

function convert(::Type{String}, cigarString::CIGARString)
    outString = ""
    for cigar in cigarString
        outString *= String(cigar)
    end
    return outString
end

function show(io::IO, cigarstr::CIGARString)
    write(io, convert(String, cigarstr))
end

function writemime(io::IO, ::MIME{symbol("text/plain")}, cs::CIGARString)
    show(io, cs)
end

function *(a::CIGARString, b::CIGARString)
    return [a;b]
end

function *(a::Operation, n::Int)
    return CIGAR(a, n)
end

function viewPosition(x::CIGARString, position::Int)
    nextIdx = start(x)
    currentPos = position
    while !done(x, nextIdx)
        currentCIGAR, nextIdx = next(x, nextIdx)
        currentPos -= currentCIGAR.Size
        if currentPos <= 0
            return nextIdx - 1
        end
    end
end
