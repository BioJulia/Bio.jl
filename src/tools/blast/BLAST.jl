module BLAST

export blastn,
       blastp,
       readblastXML,
       BLASTResult

using Bio.Seq,
      Bio.Align,
      LightXML

include("blastcommandline.jl")

end # module BLAST
