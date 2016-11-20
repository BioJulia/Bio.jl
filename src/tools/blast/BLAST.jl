module BLAST

export blastn,
       blastp,
       readblastXML,
       BLASTResult

using Bio.Seq,
      Bio.Align,
      EzXML

include("blastcommandline.jl")

end # module BLAST
