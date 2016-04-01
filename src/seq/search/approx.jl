# Approxiamte Search
# ==================

function approxsearch(seq::Sequence, pat::Sequence, k::Integer, start::Integer=1)
    return approxsearch(seq, pat, k, start, true)
end

function approxsearchindex(seq::Sequence, pat::Sequence, k::Integer, start::Integer=1)
    return approxsearch(seq, pat, k, start, true)[1]
end

function approxrsearch(seq::Sequence, pat::Sequence, k::Integer, start::Integer=endof(seq))
    return approxsearch(seq, pat, k, start, false)
end

function approxrsearchindex(seq::Sequence, pat::Sequence, k::Integer, start::Integer=endof(seq))
    return approxsearch(seq, pat, k, start, false)[1]
end

function approxsearch(seq, pat, k, start, forward)
    if k < 0
        throw(ArgumentError("k must be non-negative"))
    end
    start = Int(start)
    m = length(pat)
    if k ≥ m
        # delete/insert all
        return start:start-1
    end
    # TODO: BigInt is very slow, consider implementing "4.2 THE BLOCKS MODEL"
    T = m ≤ 64 ? UInt64 : m ≤ 128 ? UInt128 : BigInt
    stop, dist = approxsearch(T, seq, pat, k, start, forward)
    if stop == 0
        return 0:-1
    end
    start = alignback(seq, pat, dist, stop, forward)
    if forward
        return start:stop
    else
        return stop:start
    end
end

# Myers, Gene. "A fast bit-vector algorithm for approximate string matching
# based on dynamic programming." Journal of the ACM (JACM) 46.3 (1999): 395-415.
function approxsearch{T}(::Type{T}, seq, pat, k, start, forward)
    m = length(pat)
    n = length(seq)
    @assert T == BigInt || m ≤ sizeof(T) * 8

    # preprocess
    Σ = alphabet(typeof(seq))
    Peq = zeros(T, length(Σ))
    for i in 1:m
        y = forward ? pat[i] : pat[end-i+1]
        for x in Σ
            if Seq.iscompatible(x, y)
                Peq[UInt8(x)+1] |= one(T) << (i - 1)
            end
        end
    end

    Pv::T = (one(T) << m) - one(T)
    Mv::T = zero(T)
    dist::Int = m

    j::Int = start
    while (forward && j ≤ n) || (!forward && j ≥ 1)
        Eq = Peq[UInt8(seq[j])+1]
        Xv = Eq | Mv
        Xh = (((Eq & Pv) + Pv) $ Pv) | Eq

        Ph = Mv | ~(Xh | Pv)
        Mh = Pv & Xh
        if (Ph >> (m - 1)) & 1 != 0
            dist += 1
        elseif (Mh >> (m - 1)) & 1 != 0
            dist -= 1
        end

        if dist ≤ k
            return j, dist
        end

        Ph <<= 1
        Mh <<= 1
        Pv = Mh | ~(Xv | Ph)
        Mv = Ph & Xv
        j += ifelse(forward, 1, -1)
    end

    # not found
    return 0, -1
end

# run dynamic programming to get the starting position of the alignment
function alignback(seq, pat, dist, stop, forward)
    m = length(pat)
    n = length(seq)
    @assert dist < m

    H = Vector{Int}(m + 1)
    H[1] = 0
    for i in 1:m
        H[i+1] = i
    end

    j = jext = stop
    found = false
    while (forward && j ≥ max(stop - m - dist + 1, 1)) || (!forward && j ≤ min(stop + m + dist, n))
        y = seq[j]
        h_diag = H[1]
        for i in 1:m
            x = forward ? pat[end-i+1] : pat[i]
            h = min(
                H[i] + 1,
                H[i+1] + 1,
                h_diag + ifelse(Seq.iscompatible(x, y), 0, 1))
            h_diag = H[i+1]
            H[i+1] = h
        end
        if H[m+1] == dist
            jext = j
            found = true
        end
        j += ifelse(forward, -1, +1)
    end
    @assert found

    return jext
end

