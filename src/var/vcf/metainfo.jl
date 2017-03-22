# VCF MetaInfo
# ============

type MetaInfo
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # true iff values are indexed by keys (e.g. <ID=...>).
    dict::Bool
    # indexes
    tag::UnitRange{Int}
    val::UnitRange{Int}
    dictkey::Vector{UnitRange{Int}}
    dictval::Vector{UnitRange{Int}}
end

"""
    VCF.MetaInfo()

Create an unfilled VCF metainfo.
"""
function MetaInfo()
    return MetaInfo(UInt8[], 1:0, false, 1:0, 1:0, UnitRange{Int}[], UnitRange{Int}[])
end

"""
    VCF.MetaInfo(data::Vector{UInt8})

Create a VCF metainfo from `data` containing a VCF header line.
This function verifies the format and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new metainfo object.
"""
function MetaInfo(data::Vector{UInt8})
    return convert(MetaInfo, data)
end

function Base.convert(::Type{MetaInfo}, data::Vector{UInt8})
    metainfo = MetaInfo(data, 1:0, false, 1:0, 1:0, UnitRange{Int}[], UnitRange{Int}[])
    index!(metainfo)
    return metainfo
end

"""
    VCF.MetaInfo(str::AbstractString)

Create a VCF metainfo from `str` containing a VCF header line.
This function verifies the format and indexes fields for accessors.
"""
function MetaInfo(str::AbstractString)
    return convert(MetaInfo, str)
end

function Base.convert(::Type{MetaInfo}, str::AbstractString)
    return MetaInfo(convert(Vector{UInt8}, str))
end

function initialize!(metainfo::MetaInfo)
    metainfo.filled = 1:0
    metainfo.dict = false
    metainfo.tag = 1:0
    metainfo.val = 1:0
    empty!(metainfo.dictkey)
    empty!(metainfo.dictval)
    return metainfo
end

function datarange(metainfo::MetaInfo)
    return metainfo.filled
end

function isfilled(metainfo::MetaInfo)
    return !isempty(metainfo.filled)
end

function Base.:(==)(metainfo1::MetaInfo, metainfo2::MetaInfo)
    if isfilled(metainfo1) == isfilled(metainfo2) == true
        r1 = datarange(metainfo1)
        r2 = datarange(metainfo2)
        return length(r1) == length(r2) && memcmp(pointer(metainfo1.data, first(r1)), pointer(metainfo2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(metainfo1) == isfilled(metainfo2) == false
    end
end

function checkfilled(metainfo::MetaInfo)
    if !isfilled(metainfo)
        throw(ArgumentError("unfilled VCF metainfo"))
    end
end

function MetaInfo(base::MetaInfo; tag=nothing, value=nothing)
    checkfilled(base)
    buf = IOBuffer()
    print(buf, "##")
    if tag == nothing
        write(buf, base.data[base.tag])
    else
        print(buf, tag)
    end
    print(buf, '=')
    if value == nothing
        write(buf, base.data[base.val])
    elseif isa(value, String)
        print(buf, value)
    elseif isa(value, Associative) || isa(value, Vector)
        print(buf, '<')
        for (i, (key, val)) in enumerate(value)
            if i != 1
                print(buf, ',')
            end
            print(buf, key, '=')
            if needs_quote(val)
                print(buf, '"', escape_string(val), '"')
            else
                print(buf, val)
            end
        end
        print(buf, '>')
    end
    return MetaInfo(takebuf_array(buf))
end

function needs_quote(val::String)
    return contains(val, " ") || contains(val, ",") || contains(val, "\"") || contains(val, "\\")
end

function isequaltag(metainfo::MetaInfo, tag::AbstractString)
    checkfilled(metainfo)
    return length(metainfo.tag) == sizeof(tag) &&
           memcmp(pointer(metainfo.data, first(metainfo.tag)), pointer(tag), length(metainfo.tag)) == 0
end

function Bio.metainfotag(metainfo::MetaInfo)
    checkfilled(metainfo)
    return String(metainfo.data[metainfo.tag])
end

function Bio.metainfoval(metainfo::MetaInfo)
    checkfilled(metainfo)
    return String(metainfo.data[metainfo.val])
end

function Base.keys(metainfo::MetaInfo)
    checkfilled(metainfo)
    if !metainfo.dict
        throw(ArgumentError("not a dictionary"))
    end
    return [String(metainfo.data[r]) for r in metainfo.dictkey]
end

function Base.values(metainfo::MetaInfo)
    checkfilled(metainfo)
    if !metainfo.dict
        throw(ArgumentError("not a dictionary"))
    end
    vals = String[]
    for r in metainfo.dictval
        push!(vals, extract_metainfo_value(metainfo.data, r))
    end
    return vals
end

function Base.getindex(metainfo::MetaInfo, key::String)
    checkfilled(metainfo)
    if !metainfo.dict
        throw(ArgumentError("not a dictionary"))
    end
    for (i, r) in enumerate(metainfo.dictkey)
        n = length(r)
        if n == sizeof(key) && memcmp(pointer(metainfo.data, first(r)), pointer(key), n) == 0
            return extract_metainfo_value(metainfo.data, metainfo.dictval[i])
        end
    end
    throw(KeyError(key))
end

function extract_metainfo_value(data::Vector{UInt8}, range::UnitRange{Int})
    lo = first(range)
    hi = last(range)
    if length(range) ≥ 2 && data[lo] == data[hi] == UInt8('"')
        return unescape_string(String(data[lo+1:hi-1]))
    else
        return String(data[lo:hi])
    end
end

function Base.show(io::IO, metainfo::MetaInfo)
    print(io, summary(metainfo), ':')
    if isfilled(metainfo)
        println(io)
        println(io, "    tag: ", metainfotag(metainfo))
        print(io, "  value:")
        if metainfo.dict
            for (key, val) in zip(keys(metainfo), values(metainfo))
                print(io, ' ', key, "=\"", val, '"')
            end
        else
            print(io, ' ', metainfoval(metainfo))
        end
    else
        print(io, " <not filled>")
    end
end

function Base.write(io::IO, metainfo::MetaInfo)
    checkfilled(metainfo)
    return write(io, metainfo.data)
end
