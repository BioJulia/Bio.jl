# VCF Header
# ==========

immutable Header
    metainfo::Vector{MetaInfo}
    sampleID::Vector{String}

    function Header(metainfo::Vector{MetaInfo}, sampleID::Vector{String})
        return new(metainfo, sampleID)
    end
end

"""
    VCF.Header()

Create an empty VCF header.
"""
function Header()
    return Header(MetaInfo[], String[])
end

"""
    VCF.Header(metainfo::Vector, sampleID::Vector)

Create a VCF header with `metainfo` and `sampleID`.
"""
function Header(metainfo::Vector, sampleID::Vector)
    return Header(convert(Vector{MetaInfo}, metainfo), convert(Vector{String}, sampleID))
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
    return i > endof(header.metainfo)
end

function Base.next(header::Header, i)
    return header.metainfo[i], i + 1
end

function Base.find(header::Header, tag::AbstractString)
    return Base.filter(m -> isequaltag(m, tag), header.metainfo)
end

function Base.unshift!(header::Header, metainfo)
    unshift!(header.metainfo, convert(MetaInfo, metainfo))
    return header
end

function Base.push!(header::Header, metainfo)
    push!(header.metainfo, convert(MetaInfo, metainfo))
    return header
end

function Base.show(io::IO, header::Header)
    println(io, summary(header), ':')
    tags = Bio.metainfotag.(header.metainfo)
    println(io, "  metainfo tags: ", join(unique(tags), ' '))
      print(io, "     sample IDs: ", join(header.sampleID, ' '))
end

function Base.write(io::IO, header::Header)
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
