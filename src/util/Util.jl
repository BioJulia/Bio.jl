module Util

import Base:
    start,
    next,
    done,
    show,
    size,
    ==

using Bio.Seq:
    Sequence

export eachwindow,
    EachWindowIterator,
    missed

include("windows.jl")

end # module Util
