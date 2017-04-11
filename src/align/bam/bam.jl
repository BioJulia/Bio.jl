# BAM File Format
# ===============

module BAM

import BGZFStreams
import Bio.Align: SAM
import Bio: Bio, isfilled

include("bai.jl")
include("auxdata.jl")
include("reader.jl")
include("record.jl")
include("writer.jl")
include("overlap.jl")

end
