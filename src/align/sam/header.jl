# SAM Header
# ==========

immutable Header
    metainfo::Vector{MetaInfo}
end

"""
    SAM.Header()

Create an empty header.
"""
function Header()
    return Header(MetaInfo[])
end

function Base.copy(header::Header)
    return Header(header.metainfo)
end

function Base.eltype(::Type{Header})
    return MetaInfo
end

function Base.length(header::Header)
    return length(header.metainfo)
end

function Base.start(header::Header)
    return 1
end

function Base.done(header::Header, i)
    return i > length(header.metainfo)
end

function Base.next(header::Header, i)
    return header.metainfo[i], i + 1
end

"""
    find(header::Header, key::AbstractString)::Vector{MetaInfo}

Find metainfo objects satisfying `SAM.tag(metainfo) == key`.
"""
function Base.find(header::Header, key::AbstractString)::Vector{MetaInfo}
    return filter(m -> isequalkey(m, key), header.metainfo)
end

function Base.unshift!(header::Header, metainfo::MetaInfo)
    unshift!(header.metainfo, metainfo)
    return header
end

function Base.push!(header::Header, metainfo::MetaInfo)
    push!(header.metainfo, metainfo)
    return header
end
