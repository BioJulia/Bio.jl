# VCF MetaInfo
# ============

type VCFMetaInfo
    # data is supposed to be filled or not
    filled::Bool
    # true iff values are indexed by keys (e.g. <ID=...>).
    dict::Bool
    data::Vector{UInt8}
    key::UnitRange{Int}
    val::UnitRange{Int}
    dictkey::Vector{UnitRange{Int}}
    dictval::Vector{UnitRange{Int}}
end

function VCFMetaInfo(data::Vector{UInt8}=UInt8[])
    metainfo = VCFMetaInfo(false, false, data, 0:-1, 0:-1, [], [])
    if !isempty(data)
        index!(metainfo)
    end
    return metainfo
end

function initialize!(metainfo::VCFMetaInfo)
    metainfo.filled = false
    metainfo.dict = false
    metainfo.key = 0:-1
    metainfo.val = 0:-1
    empty!(metainfo.dictkey)
    empty!(metainfo.dictval)
    return metainfo
end

function isfilled(metainfo::VCFMetaInfo)
    return metainfo.filled
end

function checkfilled(metainfo::VCFMetaInfo)
    if !isfilled(metainfo)
        throw(ArgumentError("unfilled VCF metainfo"))
    end
end

function VCFMetaInfo(base::VCFMetaInfo; key=nothing, value=nothing)
    checkfilled(base)
    buf = IOBuffer()
    print(buf, "##")
    if key == nothing
        write(buf, base.data[base.key])
    else
        print(buf, key)
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
    return VCFMetaInfo(takebuf_array(buf))
end

function needs_quote(val::String)
    return contains(val, " ") || contains(val, ",") || contains(val, "\"") || contains(val, "\\")
end

function metainfokey(metainfo::VCFMetaInfo)
    checkfilled(metainfo)
    return String(metainfo.data[metainfo.key])
end

function metainfoval(metainfo::VCFMetaInfo)
    checkfilled(metainfo)
    return String(metainfo.data[metainfo.val])
end

function Base.keys(metainfo::VCFMetaInfo)
    checkfilled(metainfo)
    if !metainfo.dict
        throw(ArgumentError("not a dictionary"))
    end
    return [String(metainfo.data[r]) for r in metainfo.dictkey]
end

function Base.values(metainfo::VCFMetaInfo)
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

function Base.getindex(metainfo::VCFMetaInfo, key::String)
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

function Base.show(io::IO, metainfo::VCFMetaInfo)
    print(io, summary(metainfo), ':')
    if isfilled(metainfo)
        println(io)
        println(io, "    key: ", metainfokey(metainfo))
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

function Base.write(io::IO, metainfo::VCFMetaInfo)
    checkfilled(metainfo)
    return write(io, metainfo.data)
end

