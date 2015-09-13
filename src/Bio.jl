module Bio

using Compat, Docile
using Docile: @doc, @doc_str

abstract FileFormat

include("bufferedreader.jl")
include("ragel.jl")

include("seq/seq.jl")
include("intervals/intervals.jl")
include("tools/tools.jl")

end # module
