# SAM Header
# ==========

type SAMHeader
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
