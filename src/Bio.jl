module Bio

using Docile
using Docile: @doc, @doc_str

abstract FileFormat

include("bufferedreader.jl")
include("ragel.jl")

include("seq/seq.jl")
include("intervals/intervals.jl")

end # module
