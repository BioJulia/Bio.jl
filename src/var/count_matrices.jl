
type GenotypeCounts
    mat::Matrix{UInt8}
    alleles::Vector{Symbol}
    indidx::Indexer{Int}
    locidx::Indexer{UnitRange}
    popidx::Indexer{Vector{Int}}
end

function GenotypeCounts()
    return GenotypeCounts(Matrix{UInt8}(), Vector{Symbol}(), )
end
