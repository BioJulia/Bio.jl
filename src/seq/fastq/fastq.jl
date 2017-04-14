module FASTQ

import Automa
import Automa.RegExp: @re_str
import Bio: Bio, isfilled
import BufferedStreams
import BufferedStreams: BufferedInputStream

include("quality.jl")
include("record.jl")
include("reader.jl")
include("writer.jl")

end
