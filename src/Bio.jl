module Bio

abstract FileFormat

include("bufferedreader.jl")
include("ragel.jl")

include("seq/seq.jl")
include("intervals/intervals.jl")

end # module
