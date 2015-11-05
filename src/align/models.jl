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

    * `submat`: substitution matrix
    * `gap_open`: score of opening a new gap
    * `gap_extend`: score of extending a gap
"""
type AffineGapScoreModel{T} <: AbstractScoreModel{T}
    submat::AbstractSubstitutionMatrix{T}
    gap_open::T
    gap_extend::T

    function AffineGapScoreModel(submat::AbstractSubstitutionMatrix{T}, gap_open::T, gap_extend::T)
        @assert gap_open ≤ 0 "gap_open should be non-positive"
        @assert gap_extend ≤ 0 "gap_extend should be non-positive"
        return new(submat, gap_open, gap_extend)
    end
end

function AffineGapScoreModel{T}(submat::AbstractSubstitutionMatrix{T}, gap_open, gap_extend)
    return AffineGapScoreModel{T}(submat, T(gap_open), T(gap_extend))
end

function AffineGapScoreModel{T}(submat::AbstractSubstitutionMatrix{T}; gaps...)
    gaps = Dict(gaps)

    if haskey(gaps, :gap_open)
        gap_open = gaps[:gap_open]
    elseif haskey(gaps, :gap_open_penalty)
        gap_open = -gaps[:gap_open_penalty]
    else
        error("gap_open or gap_open_penalty argument should be passed")
    end

    if haskey(gaps, :gap_extend)
        gap_extend = gaps[:gap_extend]
    elseif haskey(gaps, :gap_extend_penalty)
        gap_extend = -gaps[:gap_extend_penalty]
    else
        error("gap_extend or gap_extend_penalty argument should be passed")
    end

    return AffineGapScoreModel(submat, T(gap_open), T(gap_extend))
end

function AffineGapScoreModel{T}(submat::AbstractMatrix{T}, gap_open, gap_extend)
    return AffineGapScoreModel(SubstitutionMatrix(submat), gap_open, gap_extend)
end

function AffineGapScoreModel{T}(submat::AbstractMatrix{T}; gaps...)
    return AffineGapScoreModel(SubstitutionMatrix(submat); gaps...)
end

# easy interface
function Base.call(::Type{AffineGapScoreModel}; scores...)
    scores = Dict(scores)
    match = scores[:match]
    mismatch = scores[:mismatch]
    gap_open = scores[:gap_open]
    gap_extend = scores[:gap_extend]
    match, mismatch, gap_open, gap_extend = promote(match, mismatch, gap_open, gap_extend)
    submat = DichotomousSubstitutionMatrix(match, mismatch)
    return AffineGapScoreModel(submat, gap_open, gap_extend)
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
