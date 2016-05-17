# Score Models
# ============
#
# Alignment scoring models.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


# Score Models
# ------------

"""
Supertype of score model.
"""
abstract AbstractScoreModel{T<:Real}

"""
    AffineGapScoreModel(submat, gap_open, gap_extend)
    AffineGapScoreModel(submat, gap_open=, gap_extend=)
    AffineGapScoreModel(match=, mismatch=, gap_open=, gap_extend=)

Affine gap scoring model.

This creates an affine gap scroing model object for alignment from a
substitution matrix (`submat`), a gap opening score (`gap_open`), and a gap
extending score (`gap_extend`). A consecutive gap of length `k` has a score of
`gap_open + gap_extend * k`. Note that both of the gap scores should be
non-positive.  As a shorthand of creating a dichotomous substitution matrix,
you can write as, for example,
`AffineGapScoreModel(match=5, mismatch=-3, gap_open=-2, gap_extend=-1)`.

Example
-------

    using Bio.Seq
    using Bio.Align

    # create an affine gap scoring model from a predefined substitution
    # matrix and gap opening/extending scores.
    affinegap = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)

    # run global alignment between two amino acid sequenecs
    pairalign(GlobalAlignment(), aa"IDGAAGQQL", aa"IDGATGQL", affinegap)

See also: `SubstitutionMatrix`, `pairalign`, `CostModel`
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

# handy interface
function AffineGapScoreModel(; scores...)
    scores = Dict(scores)
    match = scores[:match]
    mismatch = scores[:mismatch]
    gap_open = scores[:gap_open]
    gap_extend = scores[:gap_extend]
    match, mismatch, gap_open, gap_extend = promote(match, mismatch, gap_open, gap_extend)
    submat = DichotomousSubstitutionMatrix(match, mismatch)
    return AffineGapScoreModel(submat, gap_open, gap_extend)
end

function Base.show(io::IO, model::AffineGapScoreModel)
    println(io, summary(model), ':')
    if isa(model.submat, DichotomousSubstitutionMatrix)
        println(io, "       match = ", model.submat.match)
        println(io, "    mismatch = ", model.submat.mismatch)
    else
        print(io, "  ")
        println(io, model.submat)
    end
    println(io, "    gap_open = ", model.gap_open)
      print(io, "  gap_extend = ", model.gap_extend)
end


# Cost Models
# -----------

"""
Supertype of cost model.
"""
abstract AbstractCostModel{T}

"""
    CostModel(submat, insertion, deletion)
    CostModel(submat, insertion=, deletion=)
    CostModel(match=, mismatch=, insertion=, deletion=)

Cost model.

This creates a cost model object for alignment from substitution matrix
(`submat`), an insertion cost (`insertion`), and a deletion cost (`deletion`).
Note that both of the insertion and deletion costs should be non-negative.  As
a shorthand of creating a dichotomous substitution matrix, you can write as,
for example, `CostModel(match=0, mismatch=1, insertion=2, deletion=2)`.

Example
-------

    using Bio.Align

    # create a cost model from a substitution matrix and indel costs
    cost = CostModel(ones(128, 128) - eye(128), insertion=.5, deletion=.5)

    # run global alignment to minimize edit distance
    pairalign(EditDistance(), "intension", "execution", cost)

See also: `SubstitutionMatrix`, `pairalign`, `AffineGapScoreModel`
"""
type CostModel{T} <: AbstractCostModel{T}
    submat::AbstractSubstitutionMatrix{T}
    insertion::T
    deletion::T

    function CostModel(submat, insertion, deletion)
        @assert insertion ≥ 0 "insertion should be non-negative"
        @assert deletion ≥ 0 " deletion should be non-negative"
        return new(submat, insertion, deletion)
    end
end

function CostModel{T}(submat::AbstractSubstitutionMatrix{T}, insertion, deletion)
    return CostModel{T}(submat, insertion, deletion)
end

function CostModel{T}(submat::AbstractSubstitutionMatrix{T}; indels...)
    indels = Dict(indels)
    if haskey(indels, :insertion)
        insertion = indels[:insertion]
    else
        error("insertion should be passed")
    end
    if haskey(indels, :deletion)
        deletion = indels[:deletion]
    else
        error("deletion should be passed")
    end
    return CostModel(submat, insertion, deletion)
end

function CostModel{T}(submat::AbstractMatrix{T}, insertion, deletion)
    return CostModel(SubstitutionMatrix(submat), insertion, deletion)
end

function CostModel{T}(submat::AbstractMatrix{T}; indels...)
    return CostModel(SubstitutionMatrix(submat); indels...)
end

# handy interface
function CostModel(; costs...)
    costs = Dict(costs)
    match = costs[:match]
    mismatch = costs[:mismatch]
    insertion = costs[:insertion]
    deletion = costs[:deletion]
    match, mismatch, insertion, deletion = promote(match, mismatch, insertion, deletion)
    submat = DichotomousSubstitutionMatrix(match, mismatch)
    return CostModel(submat, insertion, deletion)
end
