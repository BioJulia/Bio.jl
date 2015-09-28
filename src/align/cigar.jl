# =========================
# CIGARs and CIGAR Strings
# =========================


# Single CIGARs
# --------------

# Currently REGEX to split 
const CIGAR_REGEX = r"\d+[\w-=]"

immutable CIGAR
    OP::Operation
    Size::Int
end


function CIGAR(op::Char, size::Int)
    return CIGAR(Operation(op), size)
end


function convert(::Type{CIGAR}, str::AbstractString)
    return CIGAR(str[end], parse(Int, str[1:end-1]))
end


function convert(::Type{AbstractString}, cigar::CIGAR)
    return string(cigar.Size, Char(cigar.OP))
end


function show(io::IO, cigar::CIGAR)
    write(io, convert(AbstractString, cigar))
end


# CIGARS or CIGAR Strings
# ------------------------

typealias CIGARS Vector{CIGAR}


function convert(::Type{CIGARS}, str::AbstractString)
    matches = matchall(CIGAR_REGEX, str)
    cigar_string = Vector{CIGAR}(length(matches))
    @inbounds for i in 1:length(matches)
        cigar_string[i] = CIGAR(matches[i])
    end
    return cigar_string
end


macro cigar_str(str)
    return CIGARS(str)
end


function convert(::Type{AbstractString}, cigar_string::CIGARS)
    out_string = ""
    for cigar in cigar_string
        out_string *= String(cigar)
    end
    return out_string
end


function show(io::IO, cigarstr::CIGARS)
    write(io, convert(AbstractString, cigarstr))
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


function view_position(x::CIGARS, position::Int)
    nextidx = start(x)
    current_pos = position
    while !done(x, nextidx)
        current_cigar, nextidx = next(x, nextidx)
        current_pos -= current_cigar.Size
        if current_pos <= 0
            return nextidx - 1
        end
    end
end

