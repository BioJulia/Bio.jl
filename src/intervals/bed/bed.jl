module BED

import Automa
import Automa.RegExp: @re_str
import Bio
import FixedPointNumbers: N0f8
import BufferedStreams
import ColorTypes

include("record.jl")
include("reader.jl")
include("writer.jl")

end
