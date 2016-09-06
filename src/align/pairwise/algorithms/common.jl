# Common Utils
# ============
#
# Common utilities shared among several algorithms (internal use only).
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# k: gap length
function affinegap_score(k, gap_open_penalty, gap_extend_penalty)
    return -(gap_open_penalty + gap_extend_penalty * k)
end


# Trace type
# ----------

# trace type for pairwise alignment
typealias Trace UInt8

# trace bitmap
const TRACE_NONE   = 0b00000
const TRACE_MATCH  = 0b00001
const TRACE_DELETE = 0b00010
const TRACE_INSERT = 0b00100
const TRACE_EXTDEL = 0b01000
const TRACE_EXTINS = 0b10000


# Traceback
# ---------

macro start_traceback()
    esc(quote
        anchor_point = (i, j)
        op = OP_INVALID
    end)
end

macro finish_traceback()
    quote
        push!(anchors, AlignmentAnchor(anchor_point, op))
        push!(anchors, AlignmentAnchor((i, j), OP_START))
        reverse!(anchors)
        pop!(anchors)  # remove OP_INVALID
    end
end

macro anchor(ex)
    esc(quote
        if op != $ex
            push!(anchors, AlignmentAnchor(anchor_point, op))
            op = $ex
            anchor_point = (i, j)
        end
        if ismatchop(op)
            i -= 1
            j -= 1
        elseif isinsertop(op)
            i -= 1
        elseif isdeleteop(op)
            j -= 1
        else
            @assert false
        end
    end)
end
