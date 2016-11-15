# ABIF Reader
# ===========

immutable ABIFDirEntry
    name::String
    number::Int32
    element_type::Int32
    element_size::Int32
    num_elements::Int32
    data_size::Int32
    data_offset::Int32
end

"""
    AbifReader(input::IO)
    Create a data reader of the ABIF file format.
    # Arguments
    * `input`: data source
"""
type AbifReader{T<:IO} <: Bio.IO.AbstractReader
    # input stream
    input::T

    # Tags
    dirs::Vector{ABIFDirEntry}
end

function AbifReader(input::IO)
    header = read_abif_header(input)
    tags   = read_abif_tags(input, header)

    return AbifReader(input, tags)
end

function Bio.IO.stream(reader::AbifReader)
    return reader.input
end

function Base.length(a::AbifReader)
    return length(a.dirs)
end

function Base.getindex(a::AbifReader, t::AbstractString)
    tag = filter((x) -> x.name == t, a.dirs)

    if isempty(tag)
        throw(KeyError(t))
    end

    return parse_data_tag(a, first(tag))
end

function Base.getindex(a::AbifReader, t::ABIFDirEntry)
    return parse_data_tag(a, t)
end

function Base.getindex(a::AbifReader, t::Array{ABIFDirEntry})
    return [parse_data_tag(a, k) for k in t]
end

function read_abif_header(input::IO)
    seekstart(input)
    signature = read(input, 4)

    if issupported(signature)
        version = ntoh(first(read(input, Int16, 1)))
        return parse_directory(input, position(input))
    end
end

function read_abif_tags(input::IO, header::ABIFDirEntry)
    tags = ABIFDirEntry[]
    for index in collect(1:header.num_elements)
        start = header.data_offset + index * header.element_size
        tag = parse_directory(input, start)
        push!(tags, tag)
    end
    return tags
end

function parse_directory(input::IO, pos::Int64)
    seek(input, pos)

    name         = String(read(input, 4))
    number       = Int(ntoh(first(read(input, Int32, 1))))
    element_type = Int(ntoh(first(read(input, UInt16, 1))))
    element_size = Int(ntoh(first(read(input, Int16, 1))))
    num_elements = Int(ntoh(first(read(input, Int32, 1))))
    data_size    = Int(ntoh(first(read(input, Int32, 1))))

    if data_size <= 4
        data_offset = position(input)
    else
        data_offset = ntoh(first(read(input, Int32, 1)))
    end

    ABIFDirEntry(name, number, element_type, element_size, num_elements, data_size, data_offset)
end

function parse_data_tag(a::AbifReader, tag::ABIFDirEntry)
    if tag.data_offset > 0
        seek(a.input, tag.data_offset)

        if tag.element_type == 2
            data = String(read(a.input, tag.data_size))
            return data
        elseif tag.element_type == 4
            data = read(a.input, Int16, tag.num_elements)
            data = convert_to_int!(data)
            return data
        elseif tag.element_type == 5
            data = read(a.input, Int32, tag.num_elements)
            data = convert_to_int!(data)
            return data
        elseif tag.element_type == 7
            data = read(a.input, Float32, tag.num_elements)
            data = convert_to_float!(data)
            return data
        elseif tag.element_type == 10
            year   = Int(ntoh(first(read(a.input, Int16, 1))))
            month  = Int(first(read(a.input, UInt8, 1)))
            day    = Int(first(read(a.input, UInt8, 1)))
            return "$year-$month-$day"
        elseif tag.element_type == 11
            hour    = Int(ntoh(first(read(a.input, 1))))
            minute  = Int(ntoh(first(read(a.input, 1))))
            second  = Int(ntoh(first(read(a.input, 1))))
            hsecond = Int(ntoh(first(read(a.input, 1))))
            return "$hour:$minute:$second:$hsecond"
        elseif tag.element_type == 13
            data = Bool(first(read(a.input, tag.data_size)))
            return data
        elseif tag.element_type == 18
            data = read(a.input, tag.data_size)
            return String(data[2:end])
        elseif tag.element_type == 19
            data = read(a.input, tag.data_size)
            return String(data[1:end-1])
        end
    end
end

function convert_to_int!(data::AbstractArray)
    result = Int[]

    for d in data
        r = Int(ntoh(d))
        push!(result, r)
    end

    if length(result) <= 1
        return first(result)
    end

    return result
end

function convert_to_float!(data::AbstractArray)
    result = Float32[]
    for d in data
        r = Float32(ntoh(d))
        push!(result, r)
    end

     if length(result) <= 1
        return first(result)
     end

    return result
end

function issupported(signature::Array{UInt8})
    if signature != b"ABIF"
        error("Invalid File Signature")
    end
    return true
end
