# Needleman-Wunsch Algorithm
# ==========================
# 
# The Needleman-Wunsch algorithm.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type NeedlemanWunsch{T<:Union{Signed,AbstractFloat}}
    trace::Matrix{Trace}
    H::Vector{T}
    E::Vector{T}

    function NeedlemanWunsch(m::Integer, n::Integer)
        trace = Matrix{Trace}(m + 1, n + 1)
        fill!(trace, 0xff)
        H = Vector{T}(m + 1)
        E = Vector{T}(m)
        return new(trace, H, E)
    end

    function NeedlemanWunsch()
        return NeedlemanWunsch{T}(m, n)
    end
end

function ensureroom!(nw::NeedlemanWunsch, m, n)
    if size(nw.trace, 1) < m + 1 || size(nw.trace, 2) < n + 1
        # TODO: resize the trace matrix if possible
        nw.trace = Matrix{Trace}(m + 1, n + 1)
        nw.H = resize!(nw.H, m + 1)
        nw.E = resize!(nw.E, m)
    end
    return nw
end

# Generate a code block that updates the current score and trace.
macro update()
    esc(quote
        e = E[i]
        g = h_diag + submat[a[i],b_j]
        h = max(e, f, g)
        h_diag = H[i+1]
        H[i+1] = h
        t = trace[i+1,j+1] | ft
        e == h && (t |= TRACE_DELETE)
        f == h && (t |= TRACE_INSERT)
        g == h && (t |= TRACE_MATCH)
        trace[i+1,j+1] = t
    end)
end

macro update_nextcol(gap_init, gap_extend)
    esc(quote
        e′ = e + $gap_extend
        e = max(e′, h + $gap_init)
        E[i] = e
        et = TRACE_NONE
        e == e′ && (et |= TRACE_EXTDEL)
        trace[i+1,j+2] = et
    end)
end

macro update_nextrow(gap_init, gap_extend)
    esc(quote
        f′ = f + $gap_extend
        f = max(f′, h + $gap_init)
        ft = TRACE_NONE
        f == f′ && (ft |= TRACE_EXTINS)
    end)
end


#    start   middle  end gaps of a
#    |       |       |
# a: ----ACCA---ATGTG---
# b: AAATACGATTGATGAGGGT

function run!{T}(
        nw::NeedlemanWunsch{T},
        a, b,
        submat::AbstractSubstitutionMatrix{T},
        # a
        start_gap_open_a::T,
        start_gap_extend_a::T,
        middle_gap_open_a::T,
        middle_gap_extend_a::T,
        end_gap_open_a::T,
        end_gap_extend_a::T,
        # b
        start_gap_open_b::T,
        start_gap_extend_b::T,
        middle_gap_open_b::T,
        middle_gap_extend_b::T,
        end_gap_open_b::T,
        end_gap_extend_b::T)

    m = length(a)
    n = length(b)
    ensureroom!(nw, m, n)
    H = nw.H
    E = nw.E
    trace = nw.trace

    start_gap_init_a = start_gap_open_a + start_gap_extend_a
    middle_gap_init_a = middle_gap_open_a + middle_gap_extend_a
    end_gap_init_a = end_gap_open_a + end_gap_extend_a

    start_gap_init_b = start_gap_open_b + start_gap_extend_b
    middle_gap_init_b = middle_gap_open_b + middle_gap_extend_b
    end_gap_init_b = end_gap_open_b + end_gap_extend_b

    # ↓: insertion
    # →: deletion
    # ↘: match/mismatch
    H[1] = T(0)
    trace[1,1] = TRACE_NONE
    if m == 0 && n == 0
        return H[1]
    elseif m == 0
        for j in 1:n
            trace[1,j+1] = TRACE_DELETE
        end
        return start_gap_open_a + start_gap_extend_a * n
    elseif n == 0
        for i in 1:m
            trace[i+1,1] = TRACE_INSERT
        end
        return start_gap_open_b + start_gap_extend_b * m
    else
        # initialize the first column of H, E, and trace
        @inbounds for i in 1:m
            H[i+1] = start_gap_open_b + start_gap_extend_b * i
            E[i] = H[i+1] + (i == m ? end_gap_init_a : middle_gap_init_a)
            trace[i+1,1] = TRACE_INSERT
            trace[i+1,2] = TRACE_NONE
        end

        # run dynamic programming column by column (except the last column)
        @inbounds for j in 1:n-1
            b_j = b[j]
            h_diag = H[1]
            H[1] = start_gap_open_a + start_gap_extend_a * j
            f = H[1] + middle_gap_init_b
            ft = TRACE_NONE
            trace[1,j+1] = TRACE_DELETE
            # inner loop along the sequence a
            for i in 1:m-1
                @update
                @update_nextcol middle_gap_init_a middle_gap_extend_a
                @update_nextrow middle_gap_init_b middle_gap_extend_b
            end
            # fill the final row
            let i = m
                @update
                @update_nextcol end_gap_init_a end_gap_extend_a
            end
        end

        # fill the last column
        let j = n
            b_j = b[j]
            h_diag = H[1]
            H[1] = start_gap_open_a + start_gap_extend_a * j
            f = H[1] + end_gap_init_b
            ft = TRACE_NONE
            trace[1,j+1] = TRACE_DELETE
            # inner loop along the sequence a
            @inbounds for i in 1:m-1
                @update
                @update_nextrow end_gap_init_b end_gap_extend_b
            end
            # finally, fill the last element of the DP matrices
            let i = m
                @update
            end
        end

        return H[m+1]
    end
end

function run!(nw::NeedlemanWunsch, a, b, submat, gap_open, gap_extend)
    return run!(nw, a, b, submat,
        # a
        gap_open,
        gap_extend,
        gap_open,
        gap_extend,
        gap_open,
        gap_extend,
        # b
        gap_open,
        gap_extend,
        gap_open,
        gap_extend,
        gap_open,
        gap_extend
    )
end

function traceback(nw::NeedlemanWunsch, a, b, endpos)
    anchors = Vector{AlignmentAnchor}()
    i, j = endpos
    was_extins = false
    was_extdel = false
    @start_traceback
    while i ≥ 1 || j ≥ 1
        t = nw.trace[i+1,j+1]
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
