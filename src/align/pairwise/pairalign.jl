# Interfaces
# ==========
#
# Interface functions for pairwise sequence alignment algorithms.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# include algorithms
include("alignment.jl")
include("result.jl")
include("algorithms/common.jl")
include("algorithms/needleman_wunsch.jl")
include("algorithms/banded_needleman_wunsch.jl")
include("algorithms/smith_waterman.jl")
include("algorithms/edit_distance.jl")
include("algorithms/hamming_distance.jl")

"""
    pairalign(type, seq, ref, model, [options...])

Run pairwise alignment between two sequences: `seq` and `ref`.

Available `type`s are:

* `GlobalAlignment()`
* `LocalAlignment()`
* `SemiGlobalAlignment()`
* `OverlapAlignment()`
* `EditDistance()`
* `LevenshteinDistance()`
* `HammingDistance()`

`GlobalAlignment`, `LocalAlignment`, `SemiGlobalAlignment`, and
`OverlapAlignment` are problem that maximizes alignment score between two
sequences.  Therefore, `model` should be an instance of `AbstractScoreModel`
(e.g. `AffineGapScoreModel`).

`EditDistance`, `LevenshteinDistance`, and `HammingDistance` are problem that
minimizes alignment cost between two sequences.  As for `EditDistance`, `model`
should be an instance of `AbstractCostModel` (e.g. `CostModel`).
`LevenshteinDistance` and `HammingDistance` have predefined a cost model,
so users cannot specify a cost model for these alignment types.

When you pass the `score_only=true` or `distance_only=true` option to
`pairalign`, the result of pairwise alignment holds alignment score/distance
only.  This may enable some algorithms to run faster than calculating full
alignment result.  Other available `options` are documented for each alignemnt
type.

Example
-------

    using Bio.Seq

    # create affine gap scoring model
    affinegap = AffineGapScoreModel(
        match=5,
        mismatch=-4,
        gap_open=-5,
        gap_extend=-3
    )

    # run global alignment between two DNA sequences
    pairalign(GlobalAlignment(), dna"AGGTAG", dna"ATTG", affinegap)

    # run local alignment between two DNA sequences
    pairalign(LocalAlignment(), dna"AGGTAG", dna"ATTG", affinegap)

    # you cannot specify a cost model in LevenshteinDistance
    pairalign(LevenshteinDistance(), dna"AGGTAG", dna"ATTG")

See also: `AffineGapScoreModel`, `CostModel`
"""
function pairalign end

function pairalign{S1,S2,T}(::GlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel{T};
                          score_only::Bool=false,
                          banded::Bool=false, lower_offset::Int=0, upper_offset::Int=0)
    m = length(a)
    n = length(b)
    if banded
        if m > n
            L = m - n + lower_offset
            U = upper_offset
        else
            L = lower_offset
            U = n - m + upper_offset
        end
        bnw = BandedNeedlemanWunsch{T}(m, n, L, U)
        score = run!(bnw, a, b, score.submat, score.gap_open, score.gap_extend)
        if score_only
            return PairwiseAlignmentResult{S1,S2}(score, true)
        else
            a′ = traceback(bnw, a, b, (m, n))
            return PairwiseAlignmentResult(score, true, a′, b)
        end
    else
        nw = NeedlemanWunsch{T}(m, n)
        score = run!(nw, a, b, score.submat, score.gap_open, score.gap_extend)
        if score_only
            return PairwiseAlignmentResult{S1,S2}(score, true)
        else
            a′ = traceback(nw, a, b, (m, n))
            return PairwiseAlignmentResult(score, true, a′, b)
        end
    end
end

function pairalign{S1,S2,T}(::SemiGlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel{T};
                            score_only::Bool=false)
    m = length(a)
    n = length(b)
    nw = NeedlemanWunsch{T}(m, n)
    gap_open = score.gap_open
    gap_extend = score.gap_extend
    score = run!(nw, a, b, score.submat,
        T(0), T(0), gap_open, gap_extend, T(0), T(0),
        gap_open, gap_extend, gap_open, gap_extend, gap_open, gap_extend,
    )
    if score_only
        return PairwiseAlignmentResult{S1,S2}(score, true)
    else
        a′ = traceback(nw, a, b, (m, n))
        return PairwiseAlignmentResult(score, true, a′, b)
    end
end

function pairalign{S1,S2,T}(::OverlapAlignment, a::S1, b::S2, score::AffineGapScoreModel{T};
                            score_only::Bool=false)
    m = length(a)
    n = length(b)
    nw = NeedlemanWunsch{T}(m, n)
    score = run!(nw, a, b, score.submat,
        T(0), T(0), score.gap_open, score.gap_extend, T(0), T(0),
        T(0), T(0), score.gap_open, score.gap_extend, T(0), T(0),
    )
    if score_only
        return PairwiseAlignmentResult{S1,S2}(score, true)
    else
        a′ = traceback(nw, a, b, (m, n))
        return PairwiseAlignmentResult(score, true, a′, b)
    end
end

function pairalign{S1,S2,T}(::LocalAlignment, a::S1, b::S2, score::AffineGapScoreModel{T};
                          score_only::Bool=false)
    sw = SmithWaterman{T}(length(a), length(b))
    score, endpos = run!(sw, a, b, score.submat, score.gap_open, score.gap_extend)
    if score_only
        return PairwiseAlignmentResult{S1,S2}(score, true)
    else
        a′ = traceback(sw, a, b, endpos)
        return PairwiseAlignmentResult(score, true, a′, b)
    end
end

function pairalign{S1,S2}(::EditDistance, a::S1, b::S2, cost::CostModel;
                          distance_only::Bool=false)
    dist, trace, endpos = edit_distance(a, b, cost.submat, cost.insertion, cost.deletion)
    if distance_only
        return PairwiseAlignmentResult{S1,S2}(dist, false)
    else
        a′ = edit_traceback(a, b, trace, endpos)
        return PairwiseAlignmentResult(dist, false, a′, b)
    end
end

function pairalign{S1,S2}(::LevenshteinDistance, a::S1, b::S2;
                          distance_only::Bool=false)
    unitcost = CostModel(
        DichotomousSubstitutionMatrix{Int}(0, 1),
        insertion=1,
        deletion=1
    )
    return pairalign(EditDistance(), a, b, unitcost, distance_only=distance_only)
end

function pairalign{S1,S2}(::HammingDistance, a::S1, b::S2;
                          distance_only::Bool=false)
    dist, anchors = hamming_distance(Int, a, b)
    if distance_only
        return PairwiseAlignmentResult{S1,S2}(dist, false)
    else
        a′ = AlignedSequence(a, anchors)
        return PairwiseAlignmentResult(dist, true, a′, b)
    end
end
