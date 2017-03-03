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

function ensure_margin!(stream::BufferedStreams.BufferedInputStream)
    if stream.position * 20 > length(stream.buffer) * 19
        BufferedStreams.shiftdata!(stream)
    end
    return nothing
end

function resize_and_copy!(dst::Vector{UInt8}, src::Vector{UInt8}, r::UnitRange{Int})
    rlen = length(r)
    if length(dst) != rlen
        resize!(dst, rlen)
    end
    copy!(dst, 1, src, first(r), rlen)
    return dst
end

end
