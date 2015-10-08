# =========================
# CIGARs and CIGAR Strings
# =========================

# A CIGAR String I think should be a lightweight intermediate between the CIGAR
# string read in from file, and the full representation of an alignment this module
# provides.

# CIGAR Strings
# -------------

const CIGAR_Regex = r"^(\d+[MN=XSHIDP])+$"

immutable CIGARString{T <: AbstractString}
    String::T
end

function CIGARString{T <: AbstractString}(str::T)
    if ismatch(CIGAR_Regex, str)
        error("Input string does not conform to CIGAR format.")
    end
    return CIGARString{T}(str)
end

# Regex string to enforce CIGAR format, and split strings.

macro cigar_str(str)
    return CIGARString(str)
end




#=
function convert(::Type{String}, cigarString::CIGARS)
    outString = ""
    for cigar in cigarString
        outString *= String(cigar)
    end
    return outString
end

function show(io::IO, cigarstr::CIGARS)
    write(io, convert(String, cigarstr))
end

function writemime(io::IO, ::MIME{symbol("text/plain")}, cs::CIGARS)
    show(io, cs)
end

function *(a::CIGARS, b::CIGARS)
    return [a;b]
end

function *(a::Operation, n::Int)
    return CIGAR(a, n)
end

function viewPosition(x::CIGARS, position::Int)
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
=#
