# SAM File Format
# ===============

module SAM

import Automa
import Automa.RegExp: @re_str
import BufferedStreams
import Bio.Exceptions: missingerror
importall Bio

include("flags.jl")
include("metainfo.jl")
include("record.jl")
include("header.jl")
include("reader.jl")
include("writer.jl")

end
