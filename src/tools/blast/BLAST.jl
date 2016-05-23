module BLAST

export blastn,
       blastp,
       readblastXML,
       BLASTResult

using Bio.Seq,
      Bio.Align,
      LightXML,
      Compat

include("blastcommandline.jl")

end # module BLAST
