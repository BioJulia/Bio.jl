module Bio

using Compat, Docile

abstract FileFormat

include("bufferedreader.jl")
include("ragel.jl")

include("seq/seq.jl")
include("intervals/intervals.jl")

end # module
