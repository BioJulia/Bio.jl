module Align

using Compat
using Base.Intrinsics
using Base.Order: Ordering

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof, ==, !=, <, >, <=, >=

import Base.Order: lt

export @cigar_str

include("operations.jl")
include("cigar.jl")
include("anchors.jl")

end
