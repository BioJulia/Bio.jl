# Bio.Align
# =========
#
# Module for sequence alignments.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Align

export
    Operation,
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
    OP_START,
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
    EDNAFULL,
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

using Bio.Seq
using Compat
import IntervalTrees

include("operations.jl")
include("anchors.jl")

include("types.jl")
include("submat.jl")
include("models.jl")
include("pairwise/pairalign.jl")

end
