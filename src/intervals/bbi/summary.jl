# BBI Summary
# ===========

# Supplemental Table 7.
immutable TotalSummary
    bases_covered::UInt64
    min_val::Float64
    max_val::Float64
    sum_data::Float64
    sum_squares::Float64
end

function Base.read(io::IO, ::Type{TotalSummary})
    return TotalSummary(
        read(io, UInt64),
        read(io, Float64), read(io, Float64),
        read(io, Float64), read(io, Float64))
end
