# Approxiamte Search
# ==================
#
# Approximate sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

function approxsearch(seq::Sequence, pat::Sequence, k::Integer,
                      start::Integer=1, stop::Integer=endof(seq))
    return _approxsearch(seq, pat, k, start, stop, true)
end

function approxsearchindex(seq::Sequence, pat::Sequence, k::Integer,
                           start::Integer=1, stop::Integer=endof(seq))
    return first(_approxsearch(seq, pat, k, start, stop, true))
end

function approxrsearch(seq::Sequence, pat::Sequence, k::Integer,
                       start::Integer=endof(seq), stop::Integer=1)
    return _approxsearch(seq, pat, k, start, stop, false)
end

function approxrsearchindex(seq::Sequence, pat::Sequence, k::Integer,
                            start::Integer=endof(seq), stop::Integer=1)
    return first(_approxsearch(seq, pat, k, start, stop, false))
end

function _approxsearch(seq, pat, k, start, stop, forward)
    # select a bit vector type
    # TODO: BigInt is very slow, consider implementing "4.2 THE BLOCKS MODEL"
    m = length(pat)
    T = m ≤ 64 ? UInt64 : m ≤ 128 ? UInt128 : BigInt

    if k ≥ m
        return start:start-1
    end

    matchstop, dist = _approxsearch(T, seq, pat, k, start, stop, forward)
    if matchstop == 0
        return 0:-1
    end

    matchstart = alignback(seq, pat, dist, start, matchstop, forward)
    if forward
        return matchstart:matchstop
    else
        return matchstop:matchstart
    end
end

# This returns the end index of a suffix sequence with up to `k` errors.
# More formally, when `forward = true`, it returns the minimum `j ∈ start:stop`
# such that `min_{g ∈ 1:j} δ(pat, seq[g:j]) ≤ k` where `δ(s, t)` is the edit
# distance between `s` and `t` sequences. See Myers' paper for details:
# Myers, Gene. "A fast bit-vector algorithm for approximate string matching
# based on dynamic programming." Journal of the ACM (JACM) 46.3 (1999): 395-415.
function _approxsearch{T}(::Type{T}, seq, pat, k, start, stop, forward)
    if k < 0
        throw(ArgumentError("the number of errors must be non-negative"))
    end

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
    dist = m
    j = forward ? max(start, 1) : min(start, n)

    if dist ≤ k
        return j, dist
    end

    while (forward && j ≤ min(stop, n)) || (!forward && j ≥ max(stop, 1))
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
            return j, dist  # found
        end

        Ph <<= 1
        Mh <<= 1
        Pv = Mh | ~(Xv | Ph)
        Mv = Ph & Xv
        j += ifelse(forward, +1, -1)
    end

    return 0, -1  # not found
end

# run dynamic programming to get the starting position of the alignment
function alignback(seq, pat, dist, start, matchstop, forward)
    m = length(pat)
    n = length(seq)

    H = Vector{Int}(m + 1)
    H[1] = 0
    for i in 1:m
        H[i+1] = i
    end

    j = ret = matchstop
    found = false
    while (forward && j ≥ max(start, 1)) || (!forward && j ≤ min(start, n))
        y = seq[j]
        h_diag = H[1]
        for i in 1:m
            x = forward ? pat[end-i+1] : pat[i]
            h = min(
                H[i] + 1,
                H[i+1] + 1,
                h_diag + ifelse(iscompatible(x, y), 0, 1))
            h_diag = H[i+1]
            H[i+1] = h
        end
        if H[m+1] == dist
            ret = j
            found = true
        end
        j += ifelse(forward, -1, +1)
    end
    @assert found

    return ret
end
