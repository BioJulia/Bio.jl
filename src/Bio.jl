module Bio

using Compat, Docile
using Docile: @doc, @doc_str

abstract FileFormat

include("stringfield.jl")
include("ragel.jl")

include("seq/seq.jl")
include("intervals/intervals.jl")

end # module
