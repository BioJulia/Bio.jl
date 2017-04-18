# BBI Header
# ==========

const BED_MAGIC = 0x8789F2EB
const WIG_MAGIC = 0x888FFC26

# Supplemental Table 5.
immutable Header
    magic::UInt32
    version::UInt16
    zoom_levels::UInt16
    chromosome_tree_offset::UInt64
    full_data_offset::UInt64
    full_index_offset::UInt64
    field_count::UInt16
    defined_field_count::UInt16
    auto_sql_offset::UInt64
    total_summary_offset::UInt64
    uncompress_buf_size::UInt32
    reserved::UInt64
end

function Base.read(io::IO, ::Type{Header})
    return Header(
        read(io, UInt32), read(io, UInt16), read(io, UInt16),
        read(io, UInt64), read(io, UInt64), read(io, UInt64),
        read(io, UInt16), read(io, UInt16), read(io, UInt64),
        read(io, UInt64), read(io, UInt32), read(io, UInt64))
end
