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

function isbedgraph(header::SectionHeader)
    return header.data_type == 0x01
end

function isvarstep(header::SectionHeader)
    return header.data_type == 0x02
end

function isfixedstep(header::SectionHeader)
    return header.data_type == 0x03
end

function Base.read(io::IO, ::Type{SectionHeader})
    return SectionHeader(
        read(io, UInt32), read(io, UInt32), read(io, UInt32),
        read(io, UInt32), read(io, UInt32),
        read(io, UInt8),  read(io, UInt8),  read(io, UInt16))
end
