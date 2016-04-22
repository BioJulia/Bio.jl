module Blast

using Bio.Seq,
      Bio.Align,
      LightXML

export blastn,
       blastp,
       readblastXML,
       BlastResult

include("blastcommandline.jl")

end # module Blast
