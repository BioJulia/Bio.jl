function edit_distance{T}(a, b, submat::AbstractSubstitutionMatrix{T}, insertion_cost::T, deletion_cost::T)
    m = length(a)
    n = length(b)
    trace = Matrix{Trace}(m + 1, n + 1)
    D = Vector{T}(m + 1)
    D[1] = T(0)
    for i in 1:m
        D[i+1] = i * insertion_cost
        trace[i+1,1] = TRACE_INSERT
    end
    for j in 1:n
        d_diag = D[1]
        D[1] = j * deletion_cost
        trace[1,j+1] = TRACE_DELETE
        for i in 1:m
            ins = D[i]   + insertion_cost
            del = D[i+1] + deletion_cost
            mat = d_diag + submat[a[i],b[j]]
            # find the best score and its trace
            best = min(del, ins, mat)
            t = TRACE_NONE
            del == best && (t |= TRACE_DELETE)
            ins == best && (t |= TRACE_INSERT)
            mat == best && (t |= TRACE_MATCH)
            d_diag = D[i+1]
            D[i+1] = best
            trace[i+1,j+1] = t
        end
    end
    return D[m+1], trace, (m, n)
end

function edit_traceback(a, b, trace, endpos)
    anchors = Vector{AlignmentAnchor}()
    i, j = endpos
    @start_traceback
    while i ≥ 1 || j ≥ 1
        t = trace[i+1,j+1]
        if t & TRACE_MATCH > 0
            if a[i] == b[j]
                @anchor OP_SEQ_MATCH
            else
                @anchor OP_SEQ_MISMATCH
            end
        elseif t & TRACE_DELETE > 0
            @anchor OP_DELETE
        elseif t & TRACE_INSERT > 0
            @anchor OP_INSERT
        else
            error("failed to trace back")
        end
    end
    @finish_traceback
    return AlignedSequence(a, anchors)
end

