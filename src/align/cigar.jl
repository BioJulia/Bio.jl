
immutable CIGAR
    OP::Operation
    Size::Int
end

function CIGAR(op::Character, size::Int)
    return CIGAR(, size)
end

function CIGAR(str::String)
    matches = matchall(r"(\d+)(\w)", str)
    cigarString = Vector{CIGAR}(length(matches))
    @inbounds for i in 1:length(matches)
        m = matches[i]
        cigarString[i] = CIGAR(m[end], parse(Int, m[1:end-1]))
    end
    return cigarString
end

macro cigar_str(str)
    return CIGAR(str)
end
