# ABIF Reader
# ===========


export get_tags, tagelements

immutable AbifDirEntry
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
    dirs::Vector{AbifDirEntry}
end

function AbifReader(input::IO)
    header = read_abif_header(input)
    tags   = read_abif_tags(input, header)

    return AbifReader(input, tags)
end

function Bio.IO.stream(reader::AbifReader)
    return reader.input
end

function Base.start(p::AbifReader)
    return 1
end

function Base.done(p::AbifReader, k)
    return k > length(p)
end

function Base.next(p::AbifReader, k)
    return getindex(p, p.dirs[k]), k + 1
end

function Base.length(a::AbifReader)
    return length(a.dirs)
end

function Base.checkbounds(p::AbifReader, k::Integer)
    if 1 ≤ k ≤ length(p)
        return true
    end
    throw(BoundsError(p, k))
end

function Base.getindex(a::AbifReader, k::Integer)
    checkbounds(a, k)
    return Dict([parse_data_tag(a, a.dirs[k])])
end

function Base.getindex(a::AbifReader, t::AbstractString)
    tag = filter((x) -> x.name == t, a.dirs)

    if isempty(tag)
        throw(KeyError(t))
    end

    if length(tag) > 1
        return Dict([parse_data_tag(a, b) for b in tag])
    end

    return Dict([parse_data_tag(a, first(tag))])
end

function Base.getindex(a::AbifReader, t::AbifDirEntry)
    return Dict([parse_data_tag(a, t)])
end

function Base.getindex(a::AbifReader, t::Array{AbifDirEntry})
    return Dict([parse_data_tag(a, k) for k in t])
end

"""
    get_tags(input::AbifReader)
Returns all existing tags
# Arguments
* `input`: AbifReader
"""
function get_tags(a::AbifReader)
    return [tag for tag in a.dirs]
end

"""
    get_tags(input::AbifReader, tag_name::AbstractString)
Returns all existing tags by name
# Arguments
* `input`: AbifReader
* `tag_name`: AbstractString
"""
function get_tags(a::AbifReader, t::AbstractString)
    return [tag for tag in a.dirs if isequal(tag.name, t)]
end

"""
    tagelements(stream::AbifReader, tag_name::AbstractString)
Returns the number of how many tags exists by the same name
# Arguments
* `stream`: AbifReader
* `tag_name`: AbstractString
"""
function tagelements(a::AbifReader, t::AbstractString)
    return length(get_tags(a, t))
end

# extract the header
function read_abif_header(input::IO)
    seekstart(input)
    signature = read(input, 4)

    if is_abif_signature(signature)
        version = ntoh(first(read(input, Int16, 1)))
        return parse_directory(input, position(input))
    else
        error("Invalid File Signature")
    end
end

# extract all tags in file
function read_abif_tags(input::IO, header::AbifDirEntry)
    tags = AbifDirEntry[]
    for index in 1:header.num_elements
        start = header.data_offset + index * header.element_size
        tag = parse_directory(input, start)
        push!(tags, tag)
    end
    return tags
end

# extract DirEtry
function parse_directory(input::IO, pos::Int64)
    seek(input, pos)

    name         = String(read(input, 4))
    number       = ntoh(read(input, Int32))
    element_type = ntoh(read(input, UInt16))
    element_size = ntoh(read(input, Int16))
    num_elements = ntoh(read(input, Int32))
    data_size    = ntoh(read(input, Int32))

    if data_size <= 4
        data_offset = position(input)
    else
        data_offset = ntoh(read(input, Int32))
    end

    AbifDirEntry(name, number, element_type, element_size, num_elements, data_size, data_offset)
end

# read bytes according to the element type, other values are unsupported or legacy.
function parse_data_tag(a::AbifReader, tag::AbifDirEntry)
    if tag.data_offset > 0
        seek(a.input, tag.data_offset)

        data = Tuple{String, Union{AbstractString, Integer}}

        if tag.element_type == 2
            data = String(read(a.input, tag.data_size))

        elseif tag.element_type == 4
            data = read(a.input, Int16, tag.num_elements)
            data = convert_to_int(data)

        elseif tag.element_type == 5
            data = read(a.input, Int32, tag.num_elements)
            data = convert_to_int(data)

        elseif tag.element_type == 7
            data = read(a.input, Float32, tag.num_elements)
            data = convert_to_float(data)

        elseif tag.element_type == 10
            year   = ntoh(read(a.input, Int16))
            month  = Int(read(a.input, UInt8))
            day    = Int(read(a.input, UInt8))
            data   = "$year-$month-$day"

        elseif tag.element_type == 11
            hour    = Int(ntoh(read(a.input, UInt8)))
            minute  = Int(ntoh(read(a.input, UInt8)))
            second  = Int(ntoh(read(a.input, UInt8)))
            hsecond = Int(ntoh(read(a.input, UInt8)))
            data    = "$hour:$minute:$second:$hsecond"

        elseif tag.element_type == 13
            data = Bool(first(read(a.input, tag.data_size)))

        elseif tag.element_type == 18
            data = String(read(a.input, tag.data_size)[2:end])

        elseif tag.element_type == 19
            data = String(read(a.input, tag.data_size)[1:end-1])

        else
            data = read(a.input, tag.data_size)
        end

        return format_data!(data, tag, tagelements(a, tag.name))
    end
end

# if TAG has more than one element, concatenate element number on name
function format_data!(data::Any, tag::AbifDirEntry, elements::Integer)
    if elements > 1
        return ("$(tag.name)$(tag.number)", data)
    end
    return ("$(tag.name)", data)
end

# convert a array of big-endian values
function convert_to_int(data::AbstractArray)
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

# convert a array of big-endian values
function convert_to_float(data::AbstractArray)
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

# Check if the file has a valid Signature
function is_abif_signature(signature::Array{UInt8})
    if signature != b"ABIF"
        return false
    end
    return true
end
