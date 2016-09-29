# StringFields
# ============
#
# UTF8-encoded mutable string type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module StringFields

export StringField

import BufferedStreams

"""
A simplistic mutable, utf8 encoded string.
"""
type StringField <: AbstractString
    data::Vector{UInt8}

    # position of used data in `data` in bytes (not characters)
    part::UnitRange{Int}
end

function StringField()
    return StringField(UInt8[], 1:0)
end

function StringField(data::Vector{UInt8})
    return StringField(data, 1:length(data))
end

function StringField(data::SubString)
    part = data.string.data[1+data.offset:data.offset+nextind(data, data.endof)-1]
    StringField(part, 1:length(data))
end

# From base unicode/utf8.jl
function Base.endof(s::StringField)
    d = s.data
    i = s.part.stop
    i == 0 && return i
    while Base.is_valid_continuation(d[i])
        i -= 1
    end
    @assert i >= s.part.start
    return i - s.part.start + 1
end

# From base unicode/utf8.jl
function Base.isvalid(s::StringField, i::Integer)
    return (1 <= i <= endof(s.data)) &&
        !Base.is_valid_continuation(s.data[s.part.start + i - 1])
end

# From base unicode/utf8.jl
function Base.next(s::StringField, i::Int)
    d = s.data
    b = d[s.part.start + i - 1]
    if Base.is_valid_continuation(b)
        throw(UnicodeError(Base.UTF_ERR_INVALID_INDEX, i, d[s.part.start + i - 1]))
    end

    trailing = Base.utf8_trailing[b+1]
    if length(s.part) < i + trailing
        return '\ufffd', i+1
    end
    c::UInt32 = 0
    for j = 1:trailing+1
        c <<= 6
        c += d[s.part.start + i - 1]
        i += 1
    end
    c -= Base.utf8_offset[trailing+1]
    return Char(c), i
end

function Base.copy!(field::StringField, data::Vector{UInt8},
               start::Integer, stop::Integer)
    n = stop - start + 1
    if length(field.data) < n
        resize!(field.data, n)
    end
    copy!(field.data, 1, data, start, n)
    field.part = 1:n
    return n
end

function Base.empty!(field::StringField)
    field.part = 1:0
    return field
end

function Base.isempty(field::StringField)
    return field.part.start > field.part.stop
end

function Base.convert(::Type{StringField}, str::String)
    return StringField(copy(str.data), 1:length(str.data))
end

function Base.convert(::Type{StringField}, str::SubString)
    return StringField(str)
end

function Base.convert(::Type{String}, field::StringField)
    return String(field.data[field.part])
end

function Base.convert(::Type{AbstractString}, field::StringField)
    return convert(String, field)
end

function Base.write(io::IO, field::StringField)
    write(io, convert(String, field))
end

function Base.show(io::IO, field::StringField)
    show(io, convert(String, field))
end

function Base.copy(field::StringField)
    data = field.data[field.part]
    return StringField(data, 1:length(data))
end

# From Base.hash over strings in hashing2.jl
function Base.hash(field::StringField, h::UInt64)
    h += Base.memhash_seed
    return ccall(Base.memhash, UInt, (Ptr{UInt8}, Csize_t, UInt32),
                 pointer(field.data, field.part.start),
                 length(field.part), h % UInt32) + h
end

function Base.:(==)(a::StringField, b::StringField)
    return length(a) == length(b) &&
        ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t),
              pointer(a.data, a.part.start), pointer(b.data, b.part.start), length(a)) == 0
end

function Base.:(==)(a::StringField, b::BufferedStreams.BufferedOutputStream)
    if a === b
        return true
    elseif length(a) == length(b)
        return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t),
                     pointer(a.data, a.part.start), b.buffer, length(a.data)) == 0
    else
        return false
    end
end

end # module StringFields
