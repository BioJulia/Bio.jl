module BLAST

export
    blastn,
    blastp,
    readblastXML,
    BLASTResult

import Bio.Align: AlignedSequence
import Bio.Seq: BioSequence, DNASequence, AminoAcidSequence
import EzXML

include("blastcommandline.jl")

end # module BLAST
