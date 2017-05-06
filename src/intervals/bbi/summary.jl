# BBI Total/Section Summary
# =========================

# Supplemental Table 7.
immutable Summary
    bases_covered::UInt64
    min_val::Float64
    max_val::Float64
    sum_data::Float64
    sum_squares::Float64
end

const SUMMARY_SIZE = 40

function Base.read(io::IO, ::Type{Summary})
    return Summary(
        read(io, UInt64),
        read(io, Float64), read(io, Float64),
        read(io, Float64), read(io, Float64))
end

function Base.write(stream::IO, summary::Summary)
    return write(
        stream,
        summary.bases_covered,
        summary.min_val,
        summary.max_val,
        summary.sum_data,
        summary.sum_squares)
end


# Summary per section.
# These are used to create the header, total summary, and data index.
immutable SectionSummary
    # data range
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32

    # file offset to data and its size
    offset::UInt64
    datasize::UInt64
end

function compute_total_sumamry(summaries::Vector{SectionSummary})
    if isempty(summaries)
        return Summary(0, NaN, NaN, 0.0, 0.0)
    end
    # NOTE: These values are the depth of coverage for bigBed.
    covered = UInt64(0)
    minval = Inf
    maxval = -Inf
    sumval = 0.0
    sumsqval = 0.0
    return Summary(covered, minval, maxval, sumval, sumsqval)
end
