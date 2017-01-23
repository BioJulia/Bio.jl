# VCF Header
# ==========

type VCFHeader
    metainfo::Vector{VCFMetaInfo}
    sampleID::Vector{String}
end

function VCFHeader()
    return VCFHeader([], [])
end

function Base.getindex(header::VCFHeader, key::AbstractString)
    return Base.filter(x -> metainfokey(x) == key, header.metainfo)
end

function Base.write(io::IO, header::VCFHeader)
    n = 0
    for metainfo in header.metainfo
        n += write(io, metainfo, '\n')
    end
    n += write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    if !isempty(header.sampleID)
        n += write(io, '\t', "FORMAT")
    end
    for id in header.sampleID
        n += write(io, '\t', id)
    end
    n += write(io, '\n')
    return n
end
