module Align

using Compat
using Base.Intrinsics
using Base.Order: Ordering, lt

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof, ==, !=, <, >, <=, >=

include("operations.jl")

end
