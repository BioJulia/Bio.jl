# GFF3 File Format
# ================

module GFF3

import Automa
import Automa.RegExp: @re_str
import BufferedStreams
import Bio.Exceptions: missingerror
import URIParser
importall Bio
import Bio.Intervals: Interval, IntervalCollection

include("record.jl")
include("reader.jl")
include("writer.jl")

end
