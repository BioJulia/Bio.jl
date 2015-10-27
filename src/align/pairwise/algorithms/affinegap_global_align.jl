function affinegap_global_align{T}(
        a, b,
        submat::AbstractSubstitutionMatrix{T},
        gap_open_penalty::T,
        gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    go = gap_open_penalty
    ge = gap_extend_penalty
    goe = go + ge
    trace = Matrix{Trace}(m + 1, n + 1)
    etrace = Vector{Trace}(m)
    H = Vector{T}(m + 1)
    E = Vector{T}(m)
    # run dynamic programming column by column
    @inbounds begin
        H[1] = T(0)
        trace[1,1] = TRACE_NONE
        for i in 1:m
            H[i+1] = affinegap_score(i, go, ge)
            E[i] = H[i+1] - goe
            trace[i+1,1] = TRACE_INSERT
            etrace[i] = TRACE_NONE
        end
        for j in 1:n
            b_j = b[j]
            h_diag = H[1]
            H[1] = affinegap_score(j, go, ge)
            f = H[1] - goe
            ft = TRACE_NONE
            trace[1,j+1] = TRACE_DELETE
            for i in 1:m
                e = E[i]
                h = h_diag + submat[a[i],b_j]
                best = max(e, f, h)
                # trace
                t = etrace[i] | ft
                e == best && (t |= TRACE_DELETE)
                f == best && (t |= TRACE_INSERT)
                h == best && (t |= TRACE_MATCH)
                trace[i+1,j+1] = t
                # update
                h′ = best - goe
                e′ = e - ge
                e = max(e′, h′)
                et = TRACE_NONE
                e == e′ && (et |= TRACE_EXTDEL)
                f′ = f - ge
                f = max(f′, h′)
                ft = TRACE_NONE
                f == f′ && (ft |= TRACE_EXTINS)
                E[i] = e
                h_diag = H[i+1]
                H[i+1] = best
                etrace[i] = et
            end
        end
    end
    # return the best score and the trace
    return H[end], trace
end

function affinegap_global_traceback(a, b, trace, endpos)
    anchors = Vector{AlignmentAnchor}()
    i, j = endpos
    was_extins = false
    was_extdel = false
    @start_traceback
    while i ≥ 1 || j ≥ 1
        t = trace[i+1,j+1]
        @assert !(was_extins && was_extdel)
        if was_extins
            was_extins = t & TRACE_EXTINS > 0
            @anchor OP_INSERT
        elseif was_extdel
            was_extdel = t & TRACE_EXTDEL > 0
            @anchor OP_DELETE
        elseif t & TRACE_INSERT > 0
            was_extins = t & TRACE_EXTINS > 0
            @anchor OP_INSERT
        elseif t & TRACE_DELETE > 0
            was_extdel = t & TRACE_EXTDEL > 0
            @anchor OP_DELETE
        elseif t & TRACE_MATCH > 0
            if a[i] == b[j]
                @anchor OP_SEQ_MATCH
            else
                @anchor OP_SEQ_MISMATCH
            end
        else
            error("failed to trace back")
        end
    end
    @finish_traceback
    return AlignedSequence(a, anchors)
end
