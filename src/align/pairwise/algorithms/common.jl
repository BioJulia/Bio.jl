# utils shared among algorithms

# k: gap length
function affinegap_score(k, gap_open_penalty, gap_extend_penalty)
    return -(gap_open_penalty + gap_extend_penalty * k)
end


# trace type for pairwise alignment
typealias Trace UInt8

# trace bitmap
const TRACE_NONE   = 0b0000
const TRACE_MATCH  = 0b0001
const TRACE_DELETE = 0b0010
const TRACE_INSERT = 0b0100


# utils for tracing back

macro start_traceback()
    esc(quote
        anchor_point = (i, j)
        op = OP_INVALID
    end)
end

macro finish_traceback()
    quote
        push!(anchors, AlignmentAnchor(anchor_point[1], anchor_point[2], op))
        push!(anchors, AlignmentAnchor(i, j, OP_START))
        reverse!(anchors)
        pop!(anchors)  # remove OP_INVALID
    end
end

macro anchor(ex)
    esc(quote
        if op != $ex
            push!(anchors, AlignmentAnchor(anchor_point[1], anchor_point[2], op))
            op = $ex
            anchor_point = (i, j)
        end
    end)
end

macro match()
    esc(quote
        @anchor OP_SEQ_MATCH
        i -= 1
        j -= 1
        continue
    end)
end

macro mismatch()
    esc(quote
        @anchor OP_SEQ_MISMATCH
        i -= 1
        j -= 1
        continue
    end)
end

macro insert()
    esc(quote
        @anchor OP_INSERT
        i -= 1
        continue
    end)
end

macro delete()
    esc(quote
        @anchor OP_DELETE
        j -= 1
        continue
    end)
end
