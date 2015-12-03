# Pairwise-Alignment Result
# -------------------------

"""
Result of pairwise alignment
"""
type PairwiseAlignmentResult{T,S1,S2}
    # alignment score/distance
    value::T
    isscore::Bool
    aln::Nullable{PairwiseAlignment{S1,S2}}
end

function PairwiseAlignmentResult(value, isscore, seq, ref)
    return PairwiseAlignmentResult(value, isscore, Nullable(PairwiseAlignment(seq, ref)))
end

function Base.call{S1,S2}(::Type{PairwiseAlignmentResult{S1,S2}}, value, isscore)
    return PairwiseAlignmentResult(value, isscore, Nullable{PairwiseAlignment{S1,S2}}())
end

# TODO: add useful queries

# accessors
score(aln::PairwiseAlignmentResult) = aln.value
distance(aln::PairwiseAlignmentResult) = aln.value
alignment(aln::PairwiseAlignmentResult) = get(aln.aln)
hasalignment(aln::PairwiseAlignmentResult) = !isnull(aln.aln)

function Base.show{T,S1,S2}(io::IO, aln::PairwiseAlignmentResult{T,S1,S2})
    println(io, "PairwiseAlignmentResult{", T, ",", S1, ",", S2, "}:")
    if aln.isscore
        print(io, "  score: ", aln.value)
    else
        print(io, "  distance: ", aln.value)
    end
    if hasalignment(aln)
        show_pairwise_alignment(io, alignment(aln))
    end
end
