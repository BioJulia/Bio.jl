# BAM Auxiliary Data
# ==================

immutable AuxData <: Associative{String,Any}
    data::Vector{UInt8}
end

function Base.getindex(aux::AuxData, tag::AbstractString)
    checkauxtag(tag)
    return getauxvalue(aux.data, 1, UInt8(tag[1]), UInt8(tag[2]))
end

function Base.length(aux::AuxData)
    data = aux.data
    p = 1
    len = 0
    while p ≤ length(data)
        len += 1
        p = next_tag_position(data, p)
    end
    return len
end

function Base.start(aux::AuxData)
    return 1
end

function Base.done(aux::AuxData, pos)
    return pos > length(aux.data)
end

function Base.next(aux::AuxData, pos)
    data = aux.data
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

function checkauxtag(tag::AbstractString)
    if sizeof(tag) != 2
        throw(ArgumentError("tag length must be 2"))
    end
end

function getauxvalue(data::Vector{UInt8}, pos::Int, t1::UInt8, t2::UInt8)
    pos = findauxtag(data, pos, t1, t2)
    if pos == 0
        throw(KeyError(String([t1, t2])))
    end
    pos, T = loadauxtype(data, pos + 2)
    _, val = loadauxvalue(data, pos, T)
    return val
end

function loadauxtype(data::Vector{UInt8}, p::Int)
    function auxtype(b)
        return (
            b == UInt8('A') ? Char  :
            b == UInt8('c') ? Int8  :
            b == UInt8('C') ? UInt8 :
            b == UInt8('s') ? Int16 :
            b == UInt8('S') ? UInt16 :
            b == UInt8('i') ? Int32 :
            b == UInt8('I') ? UInt32 :
            b == UInt8('f') ? Float32 :
            b == UInt8('Z') ? String :
            error("invalid type tag: '$(Char(b))'"))
    end
    t = data[p]
    if t == UInt8('B')
        return p + 2, Vector{auxtype(data[p+1])}
    else
        return p + 1, auxtype(t)
    end
end

function loadauxvalue{T}(data::Vector{UInt8}, p::Int, ::Type{T})
    return p + sizeof(T), unsafe_load(Ptr{T}(pointer(data, p)))
end

function loadauxvalue(data::Vector{UInt8}, p::Int, ::Type{Char})
    return p + 1, Char(unsafe_load(pointer(data, p)))
end

function loadauxvalue{T}(data::Vector{UInt8}, p::Int, ::Type{Vector{T}})
    n = unsafe_load(Ptr{Int32}(pointer(data, p)))
    p += 4
    xs = Vector{T}(n)
    unsafe_copy!(pointer(xs), Ptr{T}(pointer(data, p)), n)
    return p + n * sizeof(T), xs
end

function loadauxvalue(data::Vector{UInt8}, p::Int, ::Type{String})
    dataptr = pointer(data, p)
    endptr = ccall(:memchr, Ptr{Void}, (Ptr{Void}, Cint, Csize_t),
                   dataptr, '\0', length(data) - p + 1)
    q::Int = p + (endptr - dataptr) - 1
    return q + 2, String(data[p:q])
end

function findauxtag(data::Vector{UInt8}, pos::Int, t1::UInt8, t2::UInt8)
    while pos ≤ length(data) && !(data[pos] == t1 && data[pos+1] == t2)
        pos = next_tag_position(data, pos)
    end
    if pos > length(data)
        return 0
    else
        return pos
    end
end

# Find the starting position of a next tag in `data` after `p`.
# `(data[p], data[p+1])` is supposed to be a current tag.
function next_tag_position(data::Vector{UInt8}, p::Int)
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
                 error("invalid type tag: '$(Char(eltyp))'")
        p += 1
        n = unsafe_load(Ptr{Int32}(pointer(data, p)))
        p += 4 + elsize * n
    else
        error("invalid type tag: '$(Char(typ))'")
    end
    return p
end
