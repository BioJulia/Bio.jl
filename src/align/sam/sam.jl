# SAM File Format
# ===============

module SAM

import Automa
import Automa.RegExp: @re_str
import Bio: Bio, isfilled
import Bio.Exceptions: missingerror
import Bio.RecordHelper: unsafe_parse_decimal
import BufferedStreams

include("flags.jl")
include("metainfo.jl")
include("record.jl")
include("header.jl")
include("reader.jl")
include("writer.jl")

end
