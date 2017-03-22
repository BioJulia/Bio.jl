# VCF File Format
# ===============

module VCF

import Automa
import Automa.RegExp: @re_str
import Bio: Bio, isfilled
import Bio.Exceptions: missingerror
import BufferedStreams

include("record.jl")
include("metainfo.jl")
include("header.jl")
include("reader.jl")
include("writer.jl")

end
