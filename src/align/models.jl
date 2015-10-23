# Score Models
# ------------

"""
Supertype of score model.

Every score model is a problem of finding the maximum score and its alignment.
"""
abstract AbstractScoreModel{T<:Real}

"""
Affine gap scoring model.

The gap penalty of length `k` is `gap_open_penalty + gap_extend_penalty * k`.

Fields:

    * `submat`: a substitution matrix
    * `gap_open_penalty`: a penalty of opening a new gap
    * `gap_extend_penalty`: a penalty of extending a gap
"""
type AffineGapScoreModel{T} <: AbstractScoreModel{T}
    submat::AbstractSubstitutionMatrix{T}
    gap_open_penalty::T
    gap_extend_penalty::T

    function AffineGapScoreModel(submat::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
        @assert gap_open_penalty ≥ 0 "gap_open_penalty should be non-negative"
        @assert gap_extend_penalty ≥ 0 "gap_extend_penalty should be non-negative"
        return new(submat, gap_open_penalty, gap_extend_penalty)
    end
end

# default scores (TODO: rationale)
const default_match_score = 1
const default_mismatch_score = -3
const default_gap_open_penalty = 4
const default_gap_extend_penalty = 2

function AffineGapScoreModel{T}(submat::AbstractSubstitutionMatrix{T},
                                gap_open_penalty=default_gap_open_penalty,
                                gap_extend_penalty=default_gap_extend_penalty)
    return AffineGapScoreModel{T}(submat, T(gap_open_penalty), T(gap_extend_penalty))
end

function AffineGapScoreModel{T}(submat::AbstractSubstitutionMatrix{T};
                                gap_open_penalty=default_gap_open_penalty,
                                gap_extend_penalty=default_gap_extend_penalty)
    return AffineGapScoreModel(submat, gap_open_penalty, gap_extend_penalty)
end

function AffineGapScoreModel{T}(submat::AbstractMatrix{T},
                                gap_open_penalty=default_gap_open_penalty,
                                gap_extend_penalty=default_gap_extend_penalty)
    return AffineGapScoreModel(SubstitutionMatrix(submat), gap_open_penalty, gap_extend_penalty)
end

function AffineGapScoreModel{T}(submat::AbstractMatrix{T};
                                gap_open_penalty=default_gap_open_penalty,
                                gap_extend_penalty=default_gap_extend_penalty)
    return AffineGapScoreModel(SubstitutionMatrix(submat), gap_open_penalty, gap_extend_penalty)
end

# easy interface
function Base.call(::Type{AffineGapScoreModel};
                   match=default_match_score,
                   mismatch=default_mismatch_score,
                   gap_open=-default_gap_open_penalty,
                   gap_extend=-default_gap_extend_penalty)
    match, mismatch, gap_open, gap_extend = promote(match, mismatch, gap_open, gap_extend)
    submat = DichotomousSubstitutionMatrix(match, mismatch)
    return AffineGapScoreModel(submat, -gap_open, -gap_extend)
end


# Cost Models
# -----------

"""
Supertype of cost model.

Every cost model is a problem of finding the minimum cost and its alignment.
"""
abstract AbstractCostModel{T}

"""
Cost model.

Fields:

    * `submat`: a substitution matrix
    * `insertion_cost`: a cost of inserting a character into the first sequence
    * `deletion_cost`: a cost of deleting a character from the first sequence
"""
type CostModel{T} <: AbstractCostModel{T}
    submat::AbstractSubstitutionMatrix{T}
    insertion_cost::T
    deletion_cost::T

    function CostModel(submat, insertion_cost, deletion_cost)
        @assert insertion_cost ≥ 0 "insertion_cost should be non-negative"
        @assert deletion_cost ≥ 0 " deletion_cost should be non-negative"
        return new(submat, insertion_cost, deletion_cost)
    end
end

function CostModel{T}(submat::AbstractSubstitutionMatrix{T}, insertion, deletion)
    return CostModel{T}(submat, insertion, deletion)
end

function CostModel{T}(submat::AbstractSubstitutionMatrix{T};
                      insertion=T(0), deletion=T(0))
    return CostModel(submat, insertion, deletion)
end

function CostModel{T}(submat::AbstractMatrix{T};
                      insertion=T(0), deletion=T(0))
    return CostModel(SubstitutionMatrix(submat), insertion, deletion)
end

function Base.call{T}(::Type{CostModel}; match::T=T(0), mismatch::T=T(0), insertion::T=T(0), deletion::T=T(0))
    submat = DichotomousSubstitutionMatrix(match, mismatch)
    return CostModel{T}(submat, insertion, deletion)
end
