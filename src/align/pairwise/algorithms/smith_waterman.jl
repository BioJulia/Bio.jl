# Smith-Waterman Algorithm
# ========================
#
# The Smith-Waterman algorithm.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type SmithWaterman{T<:Union{Signed,AbstractFloat}}
    trace::Matrix{Trace}
    H::Vector{T}
    E::Vector{T}

    function SmithWaterman(m::Integer, n::Integer)
        trace = Matrix{Trace}(m + 1, n + 1)
        fill!(trace, 0xff)
        H = Vector{T}(m + 1)
        E = Vector{T}(m)
        return new(trace, H, E)
    end

    function SmithWaterman()
        return SmithWaterman{T}(0, 0)
    end
end

function ensureroom!(sw::SmithWaterman, m, n)
    if size(sw.trace, 1) < m + 1 || size(sw.trace, 2) < n + 1
        # TODO: resize the trace matrix if possible
        sw.trace = Matrix{Trace}(m + 1, n + 1)
        resize!(sw.H, m + 1)
        resize!(sw.E, m)
    end
    return sw
end

function run!{T}(
        sw::SmithWaterman{T},
        a, b,
        submat::AbstractSubstitutionMatrix{T},
        gap_open_a::T,
        gap_extend_a::T,
        gap_open_b::T,
        gap_extend_b::T)

    m = length(a)
    n = length(b)
    ensureroom!(sw, m, n)
    H = sw.H
    E = sw.E
    trace = sw.trace

    gap_init_a = gap_open_a + gap_extend_a
    gap_init_b = gap_open_b + gap_extend_b

    H[1] = T(0)
    trace[1,1] = TRACE_NONE
    for i in 1:m
        H[i+1] = T(0)
        E[i] = H[i+1] + gap_init_a
        trace[i+1,1] = TRACE_NONE
        if n ≥ 1
            trace[i+1,2] = TRACE_NONE
        end
    end

    best_score = H[1]
    best_endpos = (0, 0)
    @inbounds for j in 1:n
        b_j = b[j]
        h_diag = H[1]
        H[1] = T(0)
        f = H[1] + gap_init_b
        ft = TRACE_NONE
        trace[1,j+1] = TRACE_NONE
        for i in 1:m
            e = E[i]
            g = h_diag + submat[a[i],b_j]
            h = max(T(0), e, f, g)
            h_diag = H[i+1]
            H[i+1] = h
            t = trace[i+1,j+1] | ft
            e == h && (t |= TRACE_DELETE)
            f == h && (t |= TRACE_INSERT)
            g == h && (t |= TRACE_MATCH)
            trace[i+1,j+1] = t
            if h ≥ best_score
                best_score = h
                best_endpos = (i, j)
            end
            # next E
            if j != n
                e′ = e + gap_extend_a
                e = max(e′, h + gap_init_a)
                E[i] = e
                et = TRACE_NONE
                e == e′ && (et |= TRACE_EXTDEL)
                trace[i+1,j+2] = et
            end
            # next F
            f′ = f + gap_extend_b
            f = max(f′, h + gap_init_b)
            ft = TRACE_NONE
            f == f′ && (ft |= TRACE_EXTINS)
        end
    end
    return best_score, best_endpos
end

function run!{T}(sw::SmithWaterman{T}, a, b, submat, gap_open, gap_extend)
    return run!(sw, a, b, submat, gap_open, gap_extend, gap_open, gap_extend)
end

function traceback(sw::SmithWaterman, a, b, endpos)
    anchors = Vector{AlignmentAnchor}()
    i, j = endpos
    was_extins = false
    was_extdel = false
    @start_traceback
    while sw.trace[i+1,j+1] & (TRACE_MATCH | TRACE_INSERT | TRACE_DELETE) > 0
        t = sw.trace[i+1,j+1]
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
