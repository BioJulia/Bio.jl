# Banded Needleman-Wunsch Algorithm
# =================================
#
# Banded counterpart of the Needleman-Wunsch algorithm.
#
# The banded Needleman-Wunsch algorithm is similar to the normal
# Needleman-Wunsch algorithm but search space is limited to specific cells that
# form a lower and upper bounded band in the dynamic programming matrix.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type BandedNeedlemanWunsch{T<:Union{Signed,AbstractFloat}}
    trace::Matrix{Trace}
    H::Vector{T}
    E::Vector{T}
    lower::Int
    upper::Int

    function BandedNeedlemanWunsch(
            m::Integer, n::Integer,
            lower::Integer, upper::Integer)
        lower = min(lower, m)
        upper = min(upper, n)
        width = lower + upper + 1
        trace = Matrix{Trace}(width, n + 1)
        H = Vector{T}(width)
        E = Vector{T}(width)
        return new(trace, H, E, lower, upper)
    end

    function BandedNeedlemanWunsch()
        return BandedNeedlemanWunsch{T}(0, 0, 0, 0)
    end
end

function ensureroom!(nw::BandedNeedlemanWunsch, m, n, lower, upper)
    lower = min(lower, m)
    upper = min(upper, n)
    width = lower + upper + 1
    if size(nw.trace, 1) < width || size(nw.trace, 2) < n + 1
        # TODO: resize the trace matrix if possible
        nw.trace = Matrix{Trace}(width, n + 1)
        nw.H = resize!(nw.H, width)
        nw.E = resize!(nw.E, width)
    end
    nw.lower = lower
    nw.upper = upper
    return nw
end

Base.max(x) = x

# Generate a code block that updates the current score and trace.
# For example, `@update e f g` generates:
#     # update the current score
#     e = E[i-j+U+1]
#     g = H[(i-1)-(j-1)+U+1] + submat[a[i],b_j]
#     h = max(e, f, g)
#     H[i-j1+U+1] = h
#     # update the trace
#     t = trace[i-j+U+1j+1] | ft
#     e == h && (t |= TRACE_DELETE)
#     f == h && (t |= TRACE_INSERT)
#     g == h && (t |= TRACE_MATCH)
#     trace[i-j+U+1,j+1] = t
macro update2(args...)
    @assert :g in args && length(args) ≤ 3

    ex = Expr(:block)

    # update the current score
    if :e in args
        push!(ex.args, :(e = E[i-j+U+1]))
    end
    push!(ex.args, quote
        g = H[(i-1)-(j-1)+U+1] + submat[a[i],b_j]
        h = max($(args...))
        H[i-j+U+1] = h
    end)

    # update the trace
    if :e in args && :f in args
        push!(ex.args, :(t = trace[i-j+U+1,j+1] | ft))
    elseif :e in args
        push!(ex.args, :(t = trace[i-j+U+1,j+1]))
    elseif :f in args
        push!(ex.args, :(t = ft))
    else
        push!(ex.args, :(t = TRACE_NONE))
    end
    if :e in args
        push!(ex.args, :(e == h && (t |= TRACE_DELETE)))
    end
    if :f in args
        push!(ex.args, :(f == h && (t |= TRACE_INSERT)))
    end
    push!(ex.args, quote
        g == h && (t |= TRACE_MATCH)
        trace[i-j+U+1,j+1] = t
    end)

    return esc(ex)
end

macro update2_nextcol(args...)
    if length(args) == 1
        gap_init, = args
        return esc(quote
            e = h + $gap_init
            E[i-(j+1)+U+1] = e
            et = TRACE_NONE
            trace[i-(j+1)+U+1,(j+1)+1] = et
        end)
    elseif length(args) == 2
        gap_init, gap_extend = args
        return esc(quote
            e′ = e + $gap_extend
            e = max(e′, h + $gap_init)
            E[i-(j+1)+U+1] = e
            et = TRACE_NONE
            e == e′ && (et |= TRACE_EXTDEL)
            trace[i-(j+1)+U+1,(j+1)+1] = et
        end)
    end
    @assert false
end

macro update2_nextrow(gap_init, gap_extend)
    esc(quote
        f′ = f + $gap_extend
        f = max(f′, h + $gap_init)
        ft = TRACE_NONE
        f == f′ && (ft |= TRACE_EXTINS)
    end)
end

islowerbound(i, j, L) = j == i - L

function run!{T}(
        nw::BandedNeedlemanWunsch{T},
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
    L = nw.lower
    U = nw.upper
    ensureroom!(nw, m, n, L, U)
    H = nw.H
    E = nw.E
    trace = nw.trace

    start_gap_init_a = start_gap_open_a + start_gap_extend_a
    middle_gap_init_a = middle_gap_open_a + middle_gap_extend_a
    end_gap_init_a = end_gap_open_a + end_gap_extend_a

    start_gap_init_b = start_gap_open_b + start_gap_extend_b
    middle_gap_init_b = middle_gap_open_b + middle_gap_extend_b
    end_gap_init_b = end_gap_open_b + end_gap_extend_b

    # In order to save the working space, the matrices are vertically sheared.
    # Namely, the coordinate is transformed as (i, j) → (i-j+U, j), where
    # (i, j) points to the subproblem of the optimal alignment between two
    # prefixes, a[1:i] and b[1:j], within the band.
    # The position (i, j) is accessed with:
    #   * trace[i-j+U+1,j+1]
    #   * H[i-j+U+1]
    #   * E[i-j+U+1]  (j ≥ 1)

    H[0-0+U+1] = T(0)
    trace[0-0+U+1,0+1] = TRACE_NONE
    if m == 0 && n == 0
        return H[0-0+U+1]
    elseif m == 0
        for j in 1:U
            trace[0-j+U+1,j+1] = TRACE_DELETE
        end
        return start_gap_open_a + start_gap_extend_a * n
    elseif n == 0
        for i in 1:L
            trace[i-0+U+1,0+1] = TRACE_INSERT
        end
        return start_gap_open_b + start_gap_extend_b * m
    else
        # initialize the first column of H, E, and trace
        @inbounds for i in 1:L
            H[i-0+U+1] = start_gap_open_b + start_gap_extend_b * i
            E[i-1+U+1] = H[i-0+U+1] + (i == m ? end_gap_init_a : middle_gap_init_a)
            trace[i-0+U+1,0+1] = TRACE_INSERT
            trace[i-1+U+1,1+1] = TRACE_NONE
        end

        # run dynamic programming column by column (except the last column)
        @inbounds for j in 1:n-1
            b_j = b[j]
            # vertical bounds along the j-th column
            lo = max(0, j - U)
            hi = min(m, j + L)

            # fill the first cell of the column in the band
            f, ft = let i = lo
                if lo == 0
                    h = H[i-j+U+1] = start_gap_open_a + start_gap_extend_a * j
                    trace[i-j+U+1,j+1] = TRACE_DELETE
                elseif islowerbound(i, j, L)
                    @update2 g
                else
                    @update2 e g
                end
                h + middle_gap_init_b, TRACE_NONE
            end

            for i in lo+1:hi-1
                @update2 e f g
                @update2_nextcol middle_gap_init_a middle_gap_extend_a
                @update2_nextrow middle_gap_init_b middle_gap_extend_b
            end

            # fill the last cell of the column
            if hi > lo
                let i = hi
                    et = TRACE_NONE
                    if islowerbound(i, j, L)
                        @update2 f g
                        @update2_nextcol (i == m ? end_gap_init_a : middle_gap_init_a)
                    else
                        @update2 e f g
                        @update2_nextcol (i == m ? end_gap_init_a : middle_gap_init_a) (i == m ? end_gap_extend_a : middle_gap_extend_a)
                    end
                end
            end
        end

        # fill the last column
        let j = n
            b_j = b[j]
            # vertical bounds along the j-th column
            lo = max(0, j - U)
            hi = min(m, j + L)

            # fill the first cell of the column
            f, ft = let i = lo
                if lo == 0
                    h = H[i-j+U+1] = start_gap_open_a + start_gap_extend_a * j
                    trace[i-j+U+1,j+1] = TRACE_DELETE
                elseif islowerbound(i, j, L)
                    @update2 g
                else
                    @update2 e g
                end
                h + end_gap_init_b, TRACE_NONE
            end

            @inbounds for i in lo+1:hi-1
                @update2 e f g
                @update2_nextcol end_gap_init_b end_gap_extend_b
            end

            if hi > lo
                let i = hi
                    if islowerbound(i, j, L)
                        @update2 f g
                    else
                        @update2 e f g
                    end
                end
            end
        end

        return H[m-n+U+1]
    end
end

function run!(nw::BandedNeedlemanWunsch, a, b, submat, gap_open, gap_extend)
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

function traceback(nw::BandedNeedlemanWunsch, a, b, endpos)
    U = nw.upper
    anchors = Vector{AlignmentAnchor}()
    i, j = endpos
    was_extins = false
    was_extdel = false
    @start_traceback
    while i ≥ 1 || j ≥ 1
        t = nw.trace[i-j+U+1,j+1]
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
