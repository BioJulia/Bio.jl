# Record Helper
# =============
#
# Utilities to handle records in Bio.jl.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module RecordHelper

function compare_memory(p1::Ptr, p2::Ptr, len::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, len) % Int
end

function compare_memory(data1, offset1, data2, offset2, len)
    return compare_memory(pointer(data1, offset1), pointer(data2, offset2), len)
end

# r"[0-9]+" must match `data[range]`.
function unsafe_parse_decimal{T<:Unsigned}(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int})
    x = zero(T)
    @inbounds for i in range
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return x
end

# r"[-+]?[0-9]+" must match `data[range]`.
function unsafe_parse_decimal{T<:Signed}(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int})
    lo = first(range)
    if data[lo] == UInt8('-')
        sign = T(-1)
        lo += 1
    elseif data[lo] == UInt8('+')
        sign = T(+1)
        lo += 1
    else
        sign = T(+1)
    end
    x = zero(T)
    @inbounds for i in lo:last(range)
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return sign * x
end

end
