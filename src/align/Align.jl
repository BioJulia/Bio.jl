module Align

using Base.Intrinsics
using Base.Order: Ordering
using Bio.Seq

import Base: convert,
    getindex,
    show,
    length,
    start,
    next,
    done,
    copy,
    reverse,
    show,
    endof,
    ==,
    !=,
    <,
    >,
    <=,
    >=,
    *

import IntervalTrees:
    first,
    last

import Base.Order: lt
import Base.Multimedia: MIME, writemime

export Operation,
    AlignmentAnchor,
    hasop,
    Alignment,
    AlignedSequence,
    seq2ref,
    ref2seq,
    ismatchop,
    isinsertop,
    isdeleteop,
    cigar,
    OP_MATCH,
    OP_INSERT,
    OP_DELETE,
    OP_SKIP,
    OP_SOFT_CLIP,
    OP_HARD_CLIP,
    OP_PAD,
    OP_SEQ_MATCH,
    OP_SEQ_MISMATCH,
    OP_BACK,
    OP_START

export
    # alignment types
    GlobalAlignment,
    SemiGlobalAlignment,
    OverlapAlignment,
    LocalAlignment,
    EditDistance,
    HammingDistance,
    LevenshteinDistance,
    # substitution matrices
    AbstractSubstitutionMatrix,
    SubstitutionMatrix,
    DichotomousSubstitutionMatrix,
    PAM30,
    PAM70,
    PAM250,
    BLOSUM45,
    BLOSUM50,
    BLOSUM62,
    BLOSUM80,
    BLOSUM90,
    # alignment models
    AbstractScoreModel,
    AffineGapScoreModel,
    AbstractCostModel,
    CostModel,
    # pairwise alignment
    PairwiseAlignment,
    count_matches,
    count_mismatches,
    count_insertions,
    count_deletions,
    count_aligned,
    PairwiseAlignmentResult,
    pairalign,
    score,
    distance,
    alignment,
    hasalignment


include("operations.jl")
include("anchors.jl")

include("types.jl")
include("submat.jl")
include("models.jl")
include("pairwise/pairalign.jl")

end
