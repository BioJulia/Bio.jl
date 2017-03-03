# Reader Helper
# =============
#
# Utilities to generate file readers in Bio.jl.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module ReaderHelper

import BufferedStreams

@inline function anchor!(stream::BufferedStreams.BufferedInputStream, p)
    stream.anchor = p
    stream.immobilized = true
    return stream
end

@inline function upanchor!(stream::BufferedStreams.BufferedInputStream)
    @assert stream.anchor != 0 "upanchor! called with no anchor set"
    anchor = stream.anchor
    stream.anchor = 0
    stream.immobilized = false
    return anchor
end

end
