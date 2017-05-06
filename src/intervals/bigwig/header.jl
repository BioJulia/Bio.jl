# BigBed Section Header
# =====================

# Supplemental Table 13.
immutable SectionHeader
    chrom_id::UInt32
    chrom_start::UInt32
    chrom_end::UInt32
    item_step::UInt32
    item_span::UInt32
    data_type::UInt8
    reserved::UInt8
    item_count::UInt16
end

const SECTION_HEADER_SIZE = 24

function isbedgraph(datatype::UInt8)
    return datatype == 0x01
end

function isvarstep(datatype::UInt8)
    return datatype == 0x02
end

function isfixedstep(datatype::UInt8)
    return datatype == 0x03
end

function encode_datatype(datatype::Symbol)
    if datatype == :bedgraph
        return 0x01
    elseif datatype == :varstep
        return 0x02
    elseif datatype == :fixedstep
        return 0x03
    else
        throw(ArgumentError("invalid data type: $(datatype)"))
    end
end

function Base.read(io::IO, ::Type{SectionHeader})
    return SectionHeader(
        read(io, UInt32), read(io, UInt32), read(io, UInt32),
        read(io, UInt32), read(io, UInt32),
        read(io, UInt8),  read(io, UInt8),  read(io, UInt16))
end

function Base.write(stream::IO, header::SectionHeader)
    return write(
        stream,
        header.chrom_id, header.chrom_start, header.chrom_end,
        header.item_step, header.item_span, header.data_type,
        header.reserved, header.item_count)
end
