# SAM Meta-Information
# ====================

type SAMMetaInfo
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
    tag::UnitRange{Int}
    val::UnitRange{Int}
    dictkey::Vector{UnitRange{Int}}
    dictval::Vector{UnitRange{Int}}
end

function SAMMetaInfo(data::Vector{UInt8}=UInt8[])
    metainfo = SAMMetaInfo(data, 1:0, 1:0, 1:0, UnitRange{Int}[], UnitRange{Int}[])
    if !isempty(data)
        index!(metainfo)
    end
    return metainfo
end

function SAMMetaInfo(data::AbstractString)
    return SAMMetaInfo(convert(Vector{UInt8}, data))
end

function initialize!(metainfo::SAMMetaInfo)
    metainfo.filled = 1:0
    metainfo.tag = 1:0
    metainfo.val = 1:0
    empty!(metainfo.dictkey)
    empty!(metainfo.dictval)
    return metainfo
end

function isfilled(metainfo::SAMMetaInfo)
    return !isempty(metainfo.filled)
end

function datarange(metainfo::SAMMetaInfo)
    return metainfo.filled
end

function checkfilled(metainfo::SAMMetaInfo)
    if !isfilled(metainfo)
        throw(ArgumentError("unfilled SAM metainfo"))
    end
end

function SAMMetaInfo(tag::AbstractString, value)
    buf = IOBuffer()
    if tag == "CO"  # comment
        if !isa(value, AbstractString)
            throw(ArgumentError("value must be a string"))
        end
        write(buf, "@CO\t", value)
    elseif ismatch(r"[A-Z][A-Z]", tag)
        print(buf, '@', tag)
        for (key, val) in value
            print(buf, '\t', key, ':', val)
        end
    else
        throw(ArgumentError("tag must match r\"[A-Z][A-Z]\""))
    end
    return SAMMetaInfo(takebuf_array(buf))
end

function Base.:(==)(metainfo1::SAMMetaInfo, metainfo2::SAMMetaInfo)
    if isfilled(metainfo1) == isfilled(metainfo2) == true
        r1 = datarange(metainfo1)
        r2 = datarange(metainfo2)
        return length(r1) == length(r2) && memcmp(pointer(metainfo1.data, first(r1)), pointer(metainfo2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(metainfo1) == isfilled(metainfo2) == false
    end
end

function isequalkey(metainfo::SAMMetaInfo, key::AbstractString)
    if !isfilled(metainfo) || sizeof(key) != 2
        return false
    end
    k1, k2 = UInt8(key[1]), UInt8(key[2])
    return metainfo.data[metainfo.tag[1]] == k1 && metainfo.data[metainfo.tag[2]] == k2
end

function iscomment(metainfo::SAMMetaInfo)
    return isequalkey(metainfo, "CO")
end

function Bio.metainfotag(metainfo::SAMMetaInfo)
    checkfilled(metainfo)
    return String(metainfo.data[metainfo.tag])
end

function Bio.metainfoval(metainfo::SAMMetaInfo)
    checkfilled(metainfo)
    return String(metainfo.data[metainfo.val])
end

function Base.keys(metainfo::SAMMetaInfo)
    checkfilled(metainfo)
    if iscomment(metainfo)
        throw(ArgumentError("not a dictionary"))
    end
    return [String(metainfo.data[r]) for r in metainfo.dictkey]
end

function Base.values(metainfo::SAMMetaInfo)
    checkfilled(metainfo)
    if iscomment(metainfo)
        throw(ArgumentError("not a dictionary"))
    end
    return [String(metainfo.data[r]) for r in metainfo.dictval]
end

function Base.haskey(metainfo::SAMMetaInfo, key::AbstractString)
    return findkey(metainfo, key) > 0
end

function Base.getindex(metainfo::SAMMetaInfo, key::AbstractString)
    i = findkey(metainfo, key)
    if i == 0
        throw(KeyError(key))
    end
    return String(metainfo.data[metainfo.dictval[i]])
end

function findkey(metainfo::SAMMetaInfo, key::AbstractString)
    checkfilled(metainfo)
    if sizeof(key) != 2
        return 0
    end
    t1, t2 = UInt8(key[1]), UInt8(key[2])
    for (i, k) in enumerate(metainfo.dictkey)
        if metainfo.data[first(k)] == t1 && metainfo.data[first(k)+1] == t2
            return i
        end
    end
    return 0
end

function Base.show(io::IO, metainfo::SAMMetaInfo)
    print(io, summary(metainfo), ':')
    if isfilled(metainfo)
        println(io)
        println(io, "    tag: ", metainfotag(metainfo))
          print(io, "  value:")
        if !iscomment(metainfo)
            for (key, val) in zip(keys(metainfo), values(metainfo))
                print(io, ' ', key, '=', val)
            end
        else
            print(io, ' ', metainfoval(metainfo))
        end
    else
        print(io, " <not filled>")
    end
end

function Base.print(io::IO, metainfo::SAMMetaInfo)
    write(io, metainfo)
    return nothing
end

function Base.write(io::IO, metainfo::SAMMetaInfo)
    checkfilled(metainfo)
    r = datarange(metainfo)
    return unsafe_write(io, pointer(metainfo.data, first(r)), length(r))
end
