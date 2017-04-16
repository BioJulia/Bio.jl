module BigWig

import Bio: Bio, isfilled
import Bio.Intervals: BBI
import BufferedStreams
import Libz

include("header.jl")
include("reader.jl")
include("record.jl")
include("overlap.jl")

end
