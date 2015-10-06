module Align

using Base.Intrinsics
using Base.Order: Ordering

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

import Base.Order: lt

import Base.Multimedia: MIME, writemime

export Operation,
    CIGAR,
    CIGARS,
    @cigar_str,
    AlignmentAnchor,
    hasop,
    AlignmentAnchors,
    AlignedSequence,
    OP_MATCH,
    OP_INSERT,
    OP_DELETE,
    OP_SKIP,
    OP_SOFT_CLIP,
    OP_HARD_CLIP,
    OP_PAD,
    OP_SEQ_MATCH,
    OP_SEQ_MISMATCH,
    OP_BACK


include("operations.jl")
include("anchors.jl")

end
