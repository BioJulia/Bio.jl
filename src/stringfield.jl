
"""
A simplistic mutable, utf8 encoded string.
"""
type StringField <: AbstractString
    data::Vector{UInt8}
    part::UnitRange{Int}
end


function StringField()
    return StringField(UInt8[], 1:0)
end


function Base.copy!(field::StringField, data::Vector{Uint8},
                    start::Integer, stop::Integer)
    if length(field.data) < length(data)
        resize!(field.data, length(data))
    end
    n = stop - start + 1
    copy!(field.data, 1, data, start, n)
    field.part = 1:n
    return n
end


function Base.length(field::StringField)
    return length(field.part)
end


function Base.empty!(field::StringField)
    field.part = 1:0
end


function Base.isempty(field::StringField)
    return field.part.start > field.part.stop
end


function Base.convert(::Type{UTF8String}, field::StringField)
    return UTF8String(field.data[field.part])
end


function Base.convert(::Type{String}, field::StringField)
    return convert(UTF8String, field::StringField)
end


function Base.show(io::IO, field::StringField)
    show(io, convert(UTF8String, field))
end


function Base.copy(field::StringField)
    data = field.data[field.part]
    return StringField(data, 1:length(data))
end


function Base.(:(==))(a::StringField, b::StringField)
    if a === b
        return true
    elseif length(a) == length(b)
        return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t),
                     a.data, b.data, length(a)) == 0
    else
        return false
    end
end


function Base.(:(==))(a::StringField, b::BufferedStreams.BufferedOutputStream)
    if a === b
        return true
    elseif length(a) == length(b)
        return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t),
                     a.data, b.buffer, length(a)) == 0
    else
        return false
    end
end
