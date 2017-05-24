
include("distance_types.jl")

function Bio.distance{T<:UncorrectedDistance,N}(::Type{T}, seqs::Vararg{BioSequence,N})
    counts = count(count_site(T), seqs)
    return process_count(T, counts)
end
