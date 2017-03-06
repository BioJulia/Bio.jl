# VCF Header
# ==========

immutable VCFHeader
    metainfo::Vector{VCFMetaInfo}
    sampleID::Vector{String}

    function VCFHeader(metainfo::Vector{VCFMetaInfo}, sampleID::Vector{String})
        return new(metainfo, sampleID)
    end
end

"""
    VCFHeader()

Create an empty `VCFHeader` object.
"""
function VCFHeader()
    return VCFHeader(VCFMetaInfo[], String[])
end

"""
    VCFHeader(metainfo::Vector, sampleID::Vector)

Create a `VCFHeader` object with `metainfo` and `sampleID`.
"""
function VCFHeader(metainfo::Vector, sampleID::Vector)
    return VCFHeader(convert(Vector{VCFMetaInfo}, metainfo), convert(Vector{String}, sampleID))
end

function Base.eltype(::Type{VCFHeader})
    return VCFMetaInfo
end

function Base.length(header::VCFHeader)
    return length(header.metainfo)
end

function Base.start(header::VCFHeader)
    return 1
end

function Base.done(header::VCFHeader, i)
    return i > endof(header.metainfo)
end

function Base.next(header::VCFHeader, i)
    return header.metainfo[i], i + 1
end

function Base.find(header::VCFHeader, tag::AbstractString)
    return Base.filter(m -> isequaltag(m, tag), header.metainfo)
end

function Base.unshift!(header::VCFHeader, metainfo)
    unshift!(header.metainfo, convert(VCFMetaInfo, metainfo))
    return header
end

function Base.push!(header::VCFHeader, metainfo)
    push!(header.metainfo, convert(VCFMetaInfo, metainfo))
    return header
end

function Base.show(io::IO, header::VCFHeader)
    println(io, summary(header), ':')
    tags = metainfotag.(header.metainfo)
    println(io, "  metainfo tags: ", join(unique(tags), ' '))
      print(io, "     sample IDs: ", join(header.sampleID, ' '))
end

function Base.write(io::IO, header::VCFHeader)
    n = 0
    for metainfo in header.metainfo
        n += write(io, metainfo, '\n')
    end
    n += write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    if !isempty(header.sampleID)
        n += write(io, '\t', "FORMAT")
        for id in header.sampleID
            n += write(io, '\t', id)
        end
    end
    n += write(io, '\n')
    return n
end
