module BigBed

import Automa
import Automa.RegExp: @re_str
import Bio: Bio, isfilled
import Bio.Intervals: BBI, BED
import BufferedStreams
import Libz

include("reader.jl")
include("record.jl")
include("overlap.jl")

end
