# FASTA File Format
# =================

module FASTA

import Automa
import Automa.RegExp: @re_str
import Bio: Bio, isfilled
import Bio.Exceptions: missingerror
import BufferedStreams

export description,
       identifier

include("record.jl")
include("index.jl")
include("reader.jl")
include("writer.jl")

end
