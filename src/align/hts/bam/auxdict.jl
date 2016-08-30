# AuxDataDict
# ===========
#
# Auxiliary data dictionary for BAM.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# This type is not designed for very large dictionaries: time complexities in lookup
# and update operations are O(N). The memory layout is depicted as follows:
#
#   data: |tag (2 bytes)|value type (1 byte)|value (variable)|tag|...
#
# Note that tag = 0xffff indicates the value must be skipped because it is outdated.
immutable AuxDataDict <: Associative{String,Any}
    data::Vector{UInt8}
end

function AuxDataDict()
    return AuxDataDict(UInt8[])
end

function AuxDataDict{K<:AbstractString}(pairs::Pair{K}...)
    dict = AuxDataDict()
    for (tag, val) in pairs
        dict[tag] = val
    end
    return dict
end

function Base.getindex(dict::AuxDataDict, tag::AbstractString)
    checkkeytag(tag)
    return getvalue(dict.data, 1, UInt8(tag[1]), UInt8(tag[2]))
end

function Base.setindex!(dict::AuxDataDict, val, tag::AbstractString)
    checkkeytag(tag)
    setvalue!(dict.data, 1, val, UInt8(tag[1]), UInt8(tag[2]))
    return dict
end

function Base.delete!(dict::AuxDataDict, tag::AbstractString)
    checkkeytag(tag)
    deletevalue!(dict.data, 1, UInt8(tag[1]), UInt8(tag[2]))
    return dict
end

function Base.length(dict::AuxDataDict)
    data = dict.data
    p = 1
    len = 0
    while p ≤ length(data)
        if !(data[p] == data[p+1] == 0xff)
            len += 1
        end
        p = next_tag_position(data, p)
    end
    return len
end

function Base.start(dict::AuxDataDict)
    return 1
end

function Base.done(dict::AuxDataDict, pos)
    return pos > length(dict.data)
end

function Base.next(dict::AuxDataDict, pos)
    data = dict.data
    @label doit
    t1 = data[pos]
    t2 = data[pos+1]
    pos, typ = loadauxtype(data, pos + 2)
    pos, value = loadauxvalue(data, pos, typ)
    if t1 == t2 == 0xff
        @goto doit
    end
    return Pair{String,Any}(String([t1, t2]), value), pos
end


# Internals
# ---------

function getvalue(data, pos, t1, t2)
    pos = findtag(data, pos, t1, t2)
    if pos == 0
        throw(KeyError(String([t1, t2])))
    end
    pos, T = loadauxtype(data, pos + 2)
    _, val = loadauxvalue(data, pos, T)
    return val
end

function setvalue!(data, pos, val, t1, t2)
    pos = findtag(data, pos, t1, t2)
    if pos > 0
        # cancel tag
        data[pos] = data[pos+1] = 0xff
    end

    # TODO: in-place update if possible
    pos = length(data)
    if isa(val, AbstractString)
        resize!(data, length(data) + 2 + 1 + sizeof(val) + 1)
    else
        resize!(data, length(data) + 2 + 1 + sizeof(val))
    end
    data[pos+1] = t1
    data[pos+2] = t2
    storeauxtype!(data, pos + 3, typeof(val))
    storeauxvalue!(data, pos + 4, val)
    return data
end

function deletevalue!(data, pos, t1, t2)
    pos = findtag(data, pos, t1, t2)
    if pos > 0
        # cancel tag
        data[pos] = data[pos+1] = 0xff
    end
    return data
end

function findtag(data, pos, t1, t2)
    while pos ≤ length(data) && !(data[pos] == t1 && data[pos+1] == t2)
        pos = next_tag_position(data, pos)
    end
    if pos > length(data)
        # not found
        return 0
    else
        return pos
    end
end

function loadauxtype(data, p)
    t = data[p]
    if t == UInt8('B')
        return p + 2, Vector{auxtype[data[p+1]]}
    else
        return p + 1, auxtype[t]
    end
end

function loadauxvalue{T}(data, p, ::Type{T})
    return p + sizeof(T), unsafe_load(Ptr{T}(pointer(data, p)))
end

function loadauxvalue(data, p, ::Type{Char})
    return p + 1, Char(unsafe_load(pointer(data, p)))
end

function loadauxvalue{T}(data, p, ::Type{Vector{T}})
    n = unsafe_load(Ptr{Int32}(pointer(data, p)))
    p += 4
    xs = Array(T, n)
    unsafe_copy!(pointer(xs), Ptr{T}(pointer(data, p)), n)
    return p + n * sizeof(T), xs
end

function loadauxvalue(data, p, ::Type{String})
    dataptr = pointer(data, p)
    endptr = ccall(:memchr, Ptr{Void}, (Ptr{Void}, Cint, Csize_t),
                   dataptr, '\0', length(data) - p + 1)
    q::Int = p + (endptr - dataptr) - 1
    return q + 2, String(data[p:q])
end

function storeauxtype!(data, p, T)
    if T <: AbstractVector
        data[p] = UInt8('B')
    else
        data[p] = UInt8(auxtypechar[T])
    end
    return data
end

function storeauxvalue!(data, p, val)
    if isa(val, AbstractString)
        n = length(val)
        ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Csize_t), pointer(data, p), pointer(val), n)
        data[p + n] = 0x00
    elseif isa(val, AbstractVector)
        n = length(val)
        storeauxtype!(data, p, eltype(val))
        unsafe_store!(Ptr{UInt32}(pointer(data, p + 1)), n)
        ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Csize_t), pointer(data, p + 5), pointer(val), n)
    else
        unsafe_store!(Ptr{eltype(val)}(pointer(data, p)), val)
    end
    return data
end

# Find the starting position of a next tag in `data` after `p`.
# `(data[p], data[p+1])` is supposed to be a current tag.
function next_tag_position(data, p)
    typ = Char(data[p+2])
    p += 3
    if typ == 'A'
        p += 1
    elseif typ == 'c' || typ == 'C'
        p += 1
    elseif typ == 's' || typ == 'S'
        p += 2
    elseif typ == 'i' || typ == 'I'
        p += 4
    elseif typ == 'f'
        p += 4
    elseif typ == 'd'
        p += 8
    elseif typ == 'Z' || typ == 'H'
        while data[p] != 0x00  # NULL-terminalted string
            p += 1
        end
        p += 1
    elseif typ == 'B'
        eltyp = Char(data[p])
        elsize = eltyp == 'c' || eltyp == 'C'                 ? 1 :
                 eltyp == 's' || eltyp == 'S'                 ? 2 :
                 eltyp == 'i' || eltye == 'I' || eltyp == 'f' ? 4 :
                 error("unrecognized auxiliary type: ", eltyp)
        p += 1
        n = unsafe_load(Ptr{Int32}(pointer(data, p)))
        p += elsize * n
    else
        error("unrecognized auxiliary type: ", typ)
    end
    return p
end
