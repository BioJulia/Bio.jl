# Bio.Structure
# =============
#
# Module for structural biology.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Structure

import Bio.IO: FileFormat

include("model.jl")
include("pdb.jl")
include("spatial.jl")

end # Structure
