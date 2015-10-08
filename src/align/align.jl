module Align

using Compat
using Base.Intrinsics
using Base.Order: Ordering

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof, ==, !=, <, >, <=, >=, *

import Base.Order: lt

import Base.Multimedia: MIME, writemime

export OpKind, Operation, AlignmentAnchor, hasOp,
       AlignmentAnchors, AlignedSequence, OP_GAP, OP_N, OP_MM, OP_PAD, OP_MATCH,
       OP_SCLIP, OP_MISMATCH, OP_HCLIP, OP_DELETE, OP_INSERT, OP_INVALID

include("operations.jl")
include("anchors.jl")

# I think we should have a parametric type alignment, which allows us to more genrally
# describe an alignment as consisting of a source set of sequences, and then an object
# containing the view, which might be a cigar string, or anchor array and so on.
type AlignedSequence{S, R}
    source::S
    view::R
end
