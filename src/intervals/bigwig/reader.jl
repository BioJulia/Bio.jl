# BigWig Reader
# =============

immutable Reader <: Bio.IO.AbstractReader
    stream::IO
    header::BBI.Header
    zooms::Vector{BBI.Zoom}
    summary::BBI.Summary
    # chrom name => (ID, length)
    chroms::Dict{String,Tuple{UInt32,Int}}
    # ID => chrom name
    chrom_names::Dict{UInt32,String}
    # data index
    index::BBI.RTree
end

function Base.eltype(::Type{Reader})
    return Record
end

"""
    BigWig.Reader(input::IO)

Create a reader for bigWig file format.

Note that `input` must be seekable.
"""
function Reader(input::IO)
    # read header
    header = read(input, BBI.Header)
    if header.magic != BBI.WIG_MAGIC
        error("invalid bigWig magic")
    elseif header.version < 3
        error("not a supported version of bigWig")
    end
    # read zoom objects
    zoom_headers = Vector{BBI.ZoomHeader}(header.zoom_levels)
    read!(input, zoom_headers)
    zooms = [BBI.Zoom(input, h) for h in zoom_headers]
    sort!(zooms, by=z->z.header.reduction_level)
    # read summary, B tree, and R tree
    seek(input, header.total_summary_offset)
    summary = read(input, BBI.Summary)
    btree = BBI.BTree(input, header.chromosome_tree_offset)
    chroms = Dict(name => (id, Int(len)) for (name, id, len) in BBI.chromlist(btree))
    chrom_names = Dict(id => name for (name, (id, _)) in chroms)
    index = BBI.RTree(input, header.full_index_offset)
    return Reader(input, header, zooms, summary, chroms, chrom_names, index)
end


# Record
# ------

type Record
    chromstart::UInt32
    chromend::UInt32
    value::Float32
    header::SectionHeader
    reader::Reader

    function Record(chromstart, chromend, value)
        return new(chromstart, chromend, value)
    end
end


# Iterator
# --------

type IteratorState
    stream::BufferedStreams.BufferedInputStream
    done::Bool
    header::SectionHeader
    record::Record
    n_sections::UInt64
    current_section::UInt64
    n_records::UInt16
    current_record::UInt16
end

function Base.start(reader::Reader)
    seek(reader.stream, reader.header.full_data_offset)
    # this is defined as UInt32 in the spces but actually UInt64
    section_count = read(reader.stream, UInt64)
    # dummy header
    header = SectionHeader(0, 0, 0, 0, 0, 0, 0, 0)
    return IteratorState(Libz.ZlibInflateInputStream(reader.stream), false, header, Record(), section_count, 0, 0, 0)
end

function Base.done(reader::Reader, state)
    advance!(reader, state)
    return state.done
end

function Base.next(reader::Reader, state)
    return copy(state.record), state
end

function advance!(reader::Reader, state::IteratorState)
    # find a section that has at least one record
    while state.current_section < state.n_sections && state.current_record == state.n_records
        state.header = read(state.stream, SectionHeader)
        state.current_section += 1
        state.n_records = state.header.item_count
        state.current_record = 0
    end
    if state.current_record == state.n_records
        state.done = true
        return state
    end

    # read a new record
    _read!(reader, state, state.record)
    return state
end

function _read!(reader::Reader, state, record::Record)
    @assert !state.done
    header = state.header
    if isbedgraph(header.data_type)
        chromstart = read(state.stream, UInt32)
        chromend   = read(state.stream, UInt32)
    elseif isvarstep(header.data_type)
        chromstart = read(state.stream, UInt32)
        chromend   = chromstart + header.item_span
    elseif isfixedstep(header.data_type)
        chromstart = (state.current_record == 0 ? header.chrom_start : state.record.chromstart) + header.item_step
        chromend   = chromstart + header.item_span
    else
        throw(ArgumentError("invalid data type"))
    end
    value = read(state.stream, Float32)
    record.chromstart = chromstart
    record.chromend = chromend
    record.value = value
    record.header = header
    record.reader = reader
    state.current_record += 1
    return record
end


# Statistics
# ----------

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
    mean(reader::Reader, chrom, chromstart, chromend; usezoom=false)::Float32

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
    minimum(reader::Reader, chrom, chromstart, chromend; usezoom=false)::Float32

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
    maximum(reader::Reader, chrom, chromstart, chromend; usezoom=false)::Float32

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
    minval = +Inf32
    maxval = -Inf32
    defined = false
    for record in OverlapIterator(reader, chromid, chromstart, chromend)
        minval = min(minval, record.value)
        maxval = max(maxval, record.value)
        defined = true
    end
    return defined ? (minval, maxval) : (NaN32, NaN32)
end

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
