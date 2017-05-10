# BigWig Statistics
# =================

"""
    coverage(reader, chrom, chromstart, chromend; usezoom=false)::Int

Compute the coverage of values in `[chromstart, chromend]` of `chrom`.

If `usezoom` is `true`, this function tries to use precomputed statistics (zoom)
in the file.  This is often faster but not exact in most cases.
"""
function coverage(reader::Reader, chrom::AbstractString, chromstart::Integer, chromend::Integer; usezoom=false)::Int
    chromid = reader.chroms[chrom][1]
    chromstart -= 1
    if usezoom
        zoom = BBI.find_best_zoom(reader.zooms, UInt32(chromend - chromstart))
        if !isnull(zoom)
            return BBI.coverage(get(zoom), chromid, UInt32(chromstart), UInt32(chromend))
        end
    end
    return exact_coverage(reader, chromid, UInt32(chromstart), UInt32(chromend))
end

function exact_coverage(reader::Reader, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    cov::Int = 0
    for record in OverlapIterator(reader, chromid, chromstart, chromend)
        cov += BBI.coverage2((record.chromstart, record.chromend), (chromstart, chromend))
    end
    return cov
end

"""
    mean(reader, chrom, chromstart, chromend; usezoom=false)::Float32

Compute the mean of values in `[chromstart, chromend]` of `chrom`.

This function returns `NaN32` if there are no data in that range. See `coverage`
for the `usezoom` keyword argument.
"""
function mean(reader::Reader, chrom::AbstractString, chromstart::Integer, chromend::Integer; usezoom=false)::Float32
    chromid = reader.chroms[chrom][1]
    chromstart -= 1
    if usezoom
        zoom = BBI.find_best_zoom(reader.zooms, UInt32(chromend - chromstart))
        if !isnull(zoom)
            return BBI.mean(get(zoom), chromid, UInt32(chromstart), UInt32(chromend))
        end
    end
    return exact_mean(reader, chromid, UInt32(chromstart), UInt32(chromend))
end

function exact_mean(reader::Reader, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    # compute size-weighted mean
    sum = 0.0f0
    size = 0
    for record in OverlapIterator(reader, chromid, chromstart, chromend)
        cov = BBI.coverage2((record.chromstart, record.chromend), (chromstart, chromend))
        sum += record.value * cov
        size += cov
    end
    return sum / size
end

"""
    minimum(reader, chrom, chromstart, chromend; usezoom=false)::Float32

Compute the minimum of values in `[chromstart, chromend]` of `chrom`.

This function returns `NaN32` if there are no data in that range. See `coverage`
for the `usezoom` keyword argument.
"""
function minimum(reader::Reader, chrom::AbstractString, chromstart::Integer, chromend::Integer; usezoom=false)::Float32
    chromid = reader.chroms[chrom][1]
    chromstart -= 1
    if usezoom
        zoom = BBI.find_best_zoom(reader.zooms, UInt32(chromend - chromstart))
        if !isnull(zoom)
            return BBI.minimum(get(zoom), chromid, UInt32(chromstart), UInt32(chromend))
        end
    end
    return exact_extrema(reader, chromid, UInt32(chromstart), UInt32(chromend))[1]
end

"""
    maximum(reader, chrom, chromstart, chromend; usezoom=false)::Float32

Compute the maximum of values in `[chromstart, chromend]` of `chrom`.

This function returns `NaN32` if there are no data in that range. See `coverage`
for the `usezoom` keyword argument.
"""
function maximum(reader::Reader, chrom::AbstractString, chromstart::Integer, chromend::Integer; usezoom=false)::Float32
    chromid = reader.chroms[chrom][1]
    chromstart -= 1
    if usezoom
        zoom = BBI.find_best_zoom(reader.zooms, UInt32(chromend - chromstart))
        if !isnull(zoom)
            return BBI.maximum(get(zoom), chromid, UInt32(chromstart), UInt32(chromend))
        end
    end
    return exact_extrema(reader, chromid, UInt32(chromstart), UInt32(chromend))[2]
end

function exact_extrema(reader::Reader, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    min = +Inf32
    max = -Inf32
    defined = false
    for record in OverlapIterator(reader, chromid, chromstart, chromend)
        min = Base.min(min, record.value)
        max = Base.max(max, record.value)
        defined = true
    end
    return defined ? (min, max) : (NaN32, NaN32)
end

"""
    std(reader, chrom, chromstart, chromend; usezoom=false)::Float32

Compute the standard deviation of values in `[chromstart, chromend]` of `chrom`.

See `coverage` for the `usezoom` keyword argument.
"""
function std(reader::Reader, chrom::AbstractString, chromstart::Integer, chromend::Integer; usezoom=false)::Float32
    chromid = reader.chroms[chrom][1]
    chromstart -= 1
    if usezoom
        zoom = BBI.find_best_zoom(reader.zooms, UInt32(chromend - chromstart))
        if !isnull(zoom)
            return BBI.std(get(zoom), chromid, UInt32(chromstart), UInt32(chromend))
        end
    end
    return exact_std(reader, chromid, UInt32(chromstart), UInt32(chromend))
end

function exact_std(reader::Reader, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    sum = 0.0f0
    ssq = 0.0f0
    size = 0
    for record in OverlapIterator(reader, chromid, chromstart, chromend)
        cov = BBI.coverage2((record.chromstart, record.chromend), (chromstart, chromend))
        sum += record.value * cov
        ssq += record.value^2 * cov
        size += cov
    end
    return sqrt((ssq - sum^2 / size) / (size - 1))
end
