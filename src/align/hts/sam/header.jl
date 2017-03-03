# SAM Header
# ==========

immutable SAMHeader
    metainfo::Vector{SAMMetaInfo}
end

function SAMHeader()
    return SAMHeader(SAMMetaInfo[])
end

function Base.eltype(::Type{SAMHeader})
    return SAMMetaInfo
end

function Base.length(header::SAMHeader)
    return length(header.metainfo)
end

function Base.start(header::SAMHeader)
    return 1
end

function Base.done(header::SAMHeader, i)
    return i > length(header.metainfo)
end

function Base.next(header::SAMHeader, i)
    return header.metainfo[i], i + 1
end

function Base.find(header::SAMHeader, key::AbstractString)
    return filter(m -> isequalkey(m, key), header.metainfo)
end

function Base.unshift!(header::SAMHeader, metainfo::SAMMetaInfo)
    unshift!(header.metainfo, metainfo)
    return header
end

function Base.push!(header::SAMHeader, metainfo::SAMMetaInfo)
    push!(header.metainfo, metainfo)
    return header
end
