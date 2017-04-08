# BAM File Format
# ===============

module BAM

import BGZFStreams
import Bio.Align: SAM
import Bio: Bio, isfilled

include("bai.jl")
include("auxdata.jl")
include("record.jl")
include("reader.jl")
include("writer.jl")
include("overlap.jl")

end
