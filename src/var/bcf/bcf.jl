# BCF
# ===

module BCF

import Bio: Bio, isfilled
import BGZFStreams
import BufferedStreams

include("record.jl")
include("reader.jl")
include("writer.jl")

end
