# Interfaces
# ----------

# include algorithms
include("result.jl")
include("algorithms/common.jl")
include("algorithms/affinegap_global_align.jl")
include("algorithms/affinegap_banded_global_align.jl")
include("algorithms/affinegap_local_align.jl")


function pairalign{S1,S2}(::GlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false,
                          banded::Bool=false, lower::Int=0, upper::Int=0)
    submat = score.submat
    gop = score.gap_open_penalty
    gep = score.gap_extend_penalty
    if banded
        L = lower
        U = upper
        # check whether the starting and ending positions of the DP matrix are included in the band.
        if !isinband(0, 0, L, U, a, b)
            error("the starting position is not included in the band")
        elseif !isinband(length(a), length(b), L, U, a, b)
            error("the ending position is not included in the band")
        end
        if score_only
            score, _ = affinegap_banded_global_align(a, b, L, U, submat, gop, gep)
            return PairwiseAlignment{S1,S2}(score)
        else
            score, trace = affinegap_banded_global_align(a, b, L, U, submat, gop, gep)
            a′ = affinegap_banded_global_traceback(a, b, U, trace, (length(a), length(b)))
            return PairwiseAlignment(score, a′, b)
        end
    else
        if score_only
            score, _ = affinegap_global_align(a, b, submat, gop, gep)
            return PairwiseAlignment{S1,S2}(score)
        else
            score, trace = affinegap_global_align(a, b, submat, gop, gep)
            a′ = affinegap_global_traceback(a, b, trace, (length(a), length(b)))
            return PairwiseAlignment(score, a′, b)
        end
    end
end

function pairalign{S1,S2}(::LocalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false)
    submat = score.submat
    gop = score.gap_open_penalty
    gep = score.gap_extend_penalty
    if score_only
        score, _, _ = affinegap_local_align(a, b, submat, gop, gep)
        return PairwiseAlignment{S1,S2}(score)
    else
        score, trace, best_endpos = affinegap_local_align(a, b, submat, gop, gep)
        a′ = affine_local_traceback(a, b, trace, best_endpos)
        return PairwiseAlignment(score, a′, b)
    end
end

