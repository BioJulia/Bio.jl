abstract Orientation

immutable SequenceMajor <:Orientation end
immutable SiteMajor <: Orientation end

immutable SequenceMatrix{T,O<:Orientation}
    data::Matrix{T}
end

function SequenceMatrix{T,A<:Alphabet}(::Type{T}, ::Type{SequenceMajor}, vseq::Vector{BioSequence{A}})
    nseqs = length(vseq)
    @assert nseqs > 0 throw(ArgumentError("Vector of BioSequence{$A} is empty."))
    nsites = length(vseq[1])
    @inbounds for i in 2:nseqs
        length(vseq[i]) == nsites || throw(ArgumentError("Sequences in vseq must be of same length."))
    end
    mat = Matrix{T}(nsites, nseqs)
    @inbounds for seq in 1:nseqs, site in 1:nsites
        mat[site, seq] = convert(T, vseq[seq][site])
    end
    return SequenceMatrix{T,SequenceMajor}(mat)
end

function SequenceMatrix{T,A<:Alphabet}(::Type{T}, ::Type{SiteMajor}, vseq::Vector{BioSequence{A}})
    nseqs = length(vseq)
    @assert nseqs > 0 throw(ArgumentError("Vector of BioSequence{$A} is empty."))
    nsites = length(vseq[1])
    @inbounds for i in 2:nseqs
        length(vseq[i]) == nsites || throw(ArgumentError("Sequences in vseq must be of same length."))
    end
    mat = Matrix{T}(nseqs, nsites)
    @inbounds for seq in 1:nseqs, site in 1:nsites
        mat[seq, site] = convert(T, vseq[seq][site])
    end
    return SequenceMatrix{T,SiteMajor}(mat)
end


Base.size(A::SequenceMatrix) = size(A.data)
Base.linearindexing(::Type{SequenceMatrix}) = Base.LinearFast()
@inline function Base.getindex{N<:Nucleotide,O<:Orientation}(A::SequenceMatrix{N,O}, i::Int)
    return reinterpret(N, A.data[i])
end
@inline function Base.setindex!{N<:Nucleotide,O<:Orientation}(A::SequenceMatrix, v::N, i::Int)
    A.data[i] = reinterpret(UInt8, V)
end
