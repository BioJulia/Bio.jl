# Run dynamic programming only within a band.
# Formally, a[i] and b[j] can be matched only if -L ≤ j - i ≤ U.
function affinegap_banded_global_align{T}(
        a, b,
        L::Int, U::Int,
        submat::AbstractSubstitutionMatrix{T},
        gap_open_penalty::T,
        gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    go = gap_open_penalty
    ge = gap_extend_penalty
    goe = go + ge
    L = min(L, m)
    U = min(U, n)
    # band width
    W = L + U + 1
    trace = Matrix{Trace}(W, n + 1)
    H = Vector{T}(W)
    E = Vector{T}(W)
    # In order to save the working space, the matrices are vertically sheared.
    # Namely, the coordinate is transformed as (i, j) → (i-j+U, j), where
    # (i, j) points to the subproblem of the optimal alignment between two
    # prefixes, a[1:i] and b[1:j], within the band.
    # The position (i, j) is accessed with:
    #   * trace[i-j+U+1,j+1]
    #   * H[i-j+U+1]
    #   * E[i-j+U+1]  (j ≥ 1)
    @inbounds begin
        # (i, j) = (0, 0)
        H[0-0+U+1] = T(0)
        trace[0-0+U+1,0+1] = TRACE_NONE
        for i in 1:L
            H[i-0+U+1] = affinegap_score(i, go, ge)
            E[i-1+U+1] = H[i-0+U+1] - goe
            trace[i-0+U+1,0+1] = TRACE_INSERT
        end
        # NOTE: gap_extend_penalty is added in order to avoid overflow for integers
        minimum = typemin(T) + ge
        for j in 1:n
            b_j = b[j]
            if j ≤ U
                H[0-j+U+1] = affinegap_score(j, go, ge)
                f = H[0-j+U+1] - goe
                trace[0-j+U+1,j+1] = TRACE_DELETE
            else
                f = minimum
            end
            # vertical bounds along the j-th column
            lo = max(1, j - U)
            hi = min(m, j + L)
            for i in lo:hi
                e = ifelse(i == hi, minimum, E[i-j+U+1])
                h = H[i-j+U+1] + submat[a[i],b_j]
                best = max(e, f, h)
                # trace
                t = TRACE_NONE
                e == best && (t |= TRACE_DELETE)
                f == best && (t |= TRACE_INSERT)
                h == best && (t |= TRACE_MATCH)
                # update
                E[i-(j+1)+U+1] = max(e - ge, h - goe)
                f              = max(f - ge, h - goe)
                H[i-j+U+1] = best
                trace[i-j+U+1,j+1] = t
            end
        end
    end
    return H[m-n+U+1], trace
end

function affinegap_banded_global_traceback(a, b, U, trace, endpos)
    U = min(U, length(b))
    anchors = Vector{AlignmentAnchor}()
    i, j = endpos
    @start_traceback
    while i ≥ 1 || j ≥ 1
        t = trace[i-j+U+1,j+1]
        if t & TRACE_MATCH > 0
            if a[i] == b[j]
                @anchor OP_SEQ_MATCH
            else
                @anchor OP_SEQ_MISMATCH
            end
        elseif t & TRACE_DELETE > 0
            while trace[i-j+U+1,j+1] & TRACE_DELETE > 0
                @anchor OP_DELETE
            end
        elseif t & TRACE_INSERT > 0
            while trace[i-j+U+1,j+1] & TRACE_INSERT > 0
                @anchor OP_INSERT
            end
        else
            error("failed to trace back")
        end
    end
    @finish_traceback
    return AlignedSequence(a, anchors)
end

function isinband(i, j, L, U, a, b)
    return 0 ≤ i ≤ length(a) && 0 ≤ j ≤ length(b) && i - L ≤ j ≤ i + U
end

