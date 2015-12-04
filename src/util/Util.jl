module Util

import Base:
    start,
    next,
    done,
    show

import Bio.Seq:
    Sequence

export eachwindow,
    EachWindowIterator

include("windows.jl")

end # module Util
