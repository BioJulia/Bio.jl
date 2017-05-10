# BBI Total/Section Summary
# =========================

# Supplemental Table 7.
immutable Summary
    # base coverage
    cov::UInt64
    # minimum value
    min::Float64
    # maximum value
    max::Float64
    # sum of values
    sum::Float64
    # sum of squared values
    ssq::Float64
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
        summary.cov,
        summary.min,
        summary.max,
        summary.sum,
        summary.ssq)
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
