# FASTA File Format
# =================

module FASTA

import Automa
import Automa.RegExp: @re_str
import BufferedStreams
import Bio.Exceptions: missingerror
importall Bio

include("record.jl")
include("index.jl")
include("reader.jl")
include("writer.jl")

end
