module Align

using Compat
using Base.Intrinsics
using Base.Order: Ordering

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof, ==, !=, <, >, <=, >=, *

import Base.Order: lt

import Base.Multimedia: MIME, writemime

export Operation, @cigar_str, CIGAR, CIGARString, AlignmentAnchor, hasOp,
       AlignmentAnchors, AlignedSequence, OP_GAP, OP_N, OP_MM, OP_PAD, OP_MATCH,
       OP_SCLIP, OP_MISMATCH, OP_HCLIP, OP_DELETE, OP_INSERT, OP_INVALID

include("operations.jl")
include("cigar.jl")
include("anchors.jl")

end
