module BLAST

export blastn,
       blastp,
       readblastXML,
       BlastResult

using Bio.Seq,
      Bio.Align,
      LightXML,
      Compat

include("blastcommandline.jl")

end # module BLAST
