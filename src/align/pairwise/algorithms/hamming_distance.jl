# Hamming Distance
# ================
#
# The Hamming distance algorithm.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

function hamming_distance{T}(::Type{T}, a, b)
    m = length(a)
    @assert m == length(b)
    anchors = [AlignmentAnchor(0, 0, OP_START)]
    if m == 0
        # empty string
        return 0, anchors
    end
    was_match = a[1] == b[1]
    d = ifelse(was_match, T(0), T(1))
    for i in 2:m
        if a[i] == b[i]
            if !was_match
                push!(anchors, AlignmentAnchor(i - 1, i - 1, OP_SEQ_MISMATCH))
            end
            was_match = true
        else
            if was_match
                push!(anchors, AlignmentAnchor(i - 1, i - 1, OP_SEQ_MATCH))
            end
            d += T(1)
            was_match = false
        end
    end
    push!(anchors, AlignmentAnchor(m, m, was_match ? OP_SEQ_MATCH : OP_SEQ_MISMATCH))
    return d, anchors
end
