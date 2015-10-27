# Interfaces
# ----------

# include algorithms
include("result.jl")
include("algorithms/common.jl")
include("algorithms/needleman_wunsch.jl")
include("algorithms/banded_needleman_wunsch.jl")
include("algorithms/affinegap_global_align.jl")
include("algorithms/affinegap_banded_global_align.jl")
include("algorithms/affinegap_local_align.jl")
include("algorithms/affinegap_semiglobal_align.jl")
include("algorithms/edit_distance.jl")
include("algorithms/hamming_distance.jl")


function pairalign{S1,S2,T}(::GlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel{T};
                          score_only::Bool=false,
                          banded::Bool=false, lower_offset::Int=0, upper_offset::Int=0)
    submat = score.submat
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
        score = run!(bnw, a, b, submat, -score.gap_open_penalty, -score.gap_extend_penalty)
        if score_only
            return PairwiseAlignment{S1,S2}(score, true)
        else
            a′ = traceback(bnw, a, b, (m, n))
            return PairwiseAlignment(score, true, a′, b)
        end
    else
        nw = NeedlemanWunsch{T}(m, n)
        score = run!(nw, a, b, submat, -score.gap_open_penalty, -score.gap_extend_penalty)
        if score_only
            return PairwiseAlignment{S1,S2}(score, true)
        else
            a′ = traceback(nw, a, b, (m, n))
            return PairwiseAlignment(score, true, a′, b)
        end
    end
end

function pairalign{S1,S2,T}(::SemiGlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel{T};
                            score_only::Bool=false)
    submat = score.submat
    gap_open = -score.gap_open_penalty
    gap_extend = -score.gap_extend_penalty
    nw = NeedlemanWunsch{T}(length(a), length(b))
    score = run!(nw, a, b, score.submat,
        T(0), T(0), gap_open, gap_extend, T(0), T(0),
        gap_open, gap_extend, gap_open, gap_extend, gap_open, gap_extend,
    )
    if score_only
        return PairwiseAlignment{S1,S2}(score, true)
    else
        a′ = traceback(nw, a, b, (length(a), length(b)))
        return PairwiseAlignment(score, true, a′, b)
    end
end

function pairalign{S1,S2}(::LocalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false)
    submat = score.submat
    gop = score.gap_open_penalty
    gep = score.gap_extend_penalty
    if score_only
        score, _, _ = affinegap_local_align(a, b, submat, gop, gep)
        return PairwiseAlignment{S1,S2}(score, true)
    else
        score, trace, best_endpos = affinegap_local_align(a, b, submat, gop, gep)
        a′ = affine_local_traceback(a, b, trace, best_endpos)
        return PairwiseAlignment(score, true, a′, b)
    end
end

function pairalign{S1,S2}(::EditDistance, a::S1, b::S2, cost::CostModel;
                          distance_only::Bool=false)
    submat = cost.submat
    ins = cost.insertion_cost
    del = cost.deletion_cost
    if distance_only
        distance, _, _ = edit_distance(a, b, submat, ins, del)
        return PairwiseAlignment{S1,S2}(distance, false)
    else
        distance, trace, endpos = edit_distance(a, b, submat, ins, del)
        a′ = edit_traceback(a, b, trace, endpos)
        return PairwiseAlignment(distance, false, a′, b)
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
    if distance_only
        distance, _ = hamming_distance(Int, a, b)
        return PairwiseAlignment{S1,S2}(distance, false)
    else
        distance, anchors = hamming_distance(Int, a, b)
        a′ = AlignedSequence(a, anchors)
        return PairwiseAlignment(distance, true, a′, b)
    end
end
