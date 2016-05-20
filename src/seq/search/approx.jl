# Approxiamte Search
# ==================
#
# Approximate sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
Query type for approximate sequence search.
"""
immutable ApproximateSearchQuery{S<:Sequence}
    seq::S          # query sequence
    fPcom::Vector   # compatibility vector for forward search
    bPcom::Vector   # compatibility vector for backward search
    H::Vector{Int}  # distance vector for alignback function

    function ApproximateSearchQuery(seq::Sequence, direction::Symbol)
        if direction == :forward
            fPcom = approx_preprocess(seq, true)
            bPcom = []
        elseif direction == :backward
            fPcom = []
            bPcom = approx_preprocess(seq, false)
        elseif direction == :both
            fPcom = approx_preprocess(seq, true)
            bPcom = approx_preprocess(seq, false)
        else
            throw(ArgumentError("direction '$direction' is invalid"))
        end
        H = Vector{Int}(length(seq) + 1)
        return new(seq, fPcom, bPcom, H)
    end
end

"""
    ApproximateSearchQuery(pat::Sequence[, direction=:both])

Create an query object for approximate sequence search from the `pat` sequence.

# Arguments
* `pat`: Query sequence.
* `direction=:both`: Search direction (`:forward`, `:backward`, or `:both`).
"""
function ApproximateSearchQuery(pat::Sequence, direction::Symbol=:both)
    return ApproximateSearchQuery{typeof(pat)}(pat, direction)
end

function approx_preprocess(pat, forward)
    # select a bit vector type
    # TODO: BigInt is very slow, consider implementing "4.2 THE BLOCKS MODEL"
    m = length(pat)
    T = m ≤ 64 ? UInt64 : m ≤ 128 ? UInt128 : BigInt
    Σ = alphabet(eltype(pat))
    Pcom = zeros(T, length(Σ))
    for i in 1:m
        y = forward ? pat[i] : pat[end-i+1]
        for x in Σ
            if Seq.iscompatible(x, y)
                Pcom[UInt8(x)+1] |= one(T) << (i - 1)
            end
        end
    end
    return Pcom
end

"""
    approxsearch(seq, pat, k[, start=1[, stop=endof(seq)]])

Return the range of the first occurrence of `pat` in `seq[start:stop]` allowing
up to `k` errors; symbol comparison is done using `Bio.Seq.iscompatible`.
"""
function approxsearch(seq::Sequence, pat::Sequence, k::Integer,
                      start::Integer=1, stop::Integer=endof(seq))
    return approxsearch(seq, ApproximateSearchQuery(pat, :forward), k, start, stop)
end

function approxsearch(seq::Sequence, query::ApproximateSearchQuery, k::Integer,
                      start::Integer=1, stop::Integer=endof(seq))
    return _approxsearch(query, seq, k, start, stop, true)
end

"""
    approxrsearch(seq, pat, k[, start=endof(seq)[, stop=1]])

Return the range of the last occurrence of `pat` in `seq[stop:start]` allowing
up to `k` errors; symbol comparison is done using `Bio.Seq.iscompatible`.
"""
function approxrsearch(seq::Sequence, pat::Sequence, k::Integer,
                       start::Integer=endof(seq), stop::Integer=1)
    return approxrsearch(seq, ApproximateSearchQuery(pat, :backward), k, start, stop)
end

function approxrsearch(seq::Sequence, query::ApproximateSearchQuery, k::Integer,
                       start::Integer=endof(seq), stop::Integer=1)
    return _approxsearch(query, seq, k, start, stop, false)
end

"""
    approxsearchindex(seq, pat, k[, start=1[, stop=endof(seq)]])

Return the index of the first occurrence of `pat` in `seq[start:stop]` allowing
up to `k` errors; symbol comparison is done using `Bio.Seq.iscompatible`.
"""
function approxsearchindex(seq::Sequence, pat::Sequence, k::Integer,
                           start::Integer=1, stop::Integer=endof(seq))
    return first(approxsearch(seq, pat, k, start, stop))
end

function approxsearchindex(seq::Sequence, query::ApproximateSearchQuery, k::Integer,
                           start::Integer=1, stop::Integer=endof(seq))
    return first(approxsearch(seq, query, k, start, stop))
end

"""
    approxrsearchindex(seq, pat, k[, start=endof(seq)[, stop=1]])

Return the index of the last occurrence of `pat` in `seq[stop:start]` allowing
up to `k` errors; symbol comparison is done using `Bio.Seq.iscompatible`.
"""
function approxrsearchindex(seq::Sequence, pat::Sequence, k::Integer,
                            start::Integer=endof(seq), stop::Integer=1)
    return first(approxrsearch(seq, pat, k, start, stop))
end

function approxrsearchindex(seq::Sequence, query::ApproximateSearchQuery, k::Integer,
                            start::Integer=endof(seq), stop::Integer=1)
    return first(approxrsearch(seq, query, k, start, stop))
end

function _approxsearch(query, seq, k, start, stop, forward)
    if forward && isempty(query.fPcom)
        throw(ArgumentError("query is not preprocessed for forward search"))
    end
    if !forward && isempty(query.bPcom)
        throw(ArgumentError("query is not preprocessed for backward search"))
    end

    if k ≥ length(query.seq)
        return start:start-1
    end

    # search the approximate suffix
    matchstop, dist = search_approx_suffix(
        forward ? query.fPcom : query.bPcom,
        query.seq, seq, k, start, stop, forward)
    if matchstop == 0
        return 0:-1
    end

    # locate the starting position of the match
    matchstart = alignback!(query.H, query.seq, seq, dist, start, matchstop, forward)
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
# NOTE: `Pcom` corresponds to `Peq` in the paper.
function search_approx_suffix{T}(Pcom::Vector{T}, pat, seq, k, start, stop, forward)
    if k < 0
        throw(ArgumentError("the number of errors must be non-negative"))
    end

    m = length(pat)
    n = length(seq)
    @assert T == BigInt || m ≤ sizeof(T) * 8

    Pv::T = (one(T) << m) - one(T)
    Mv::T = zero(T)
    dist = m
    j = forward ? max(start, 1) : min(start, n)

    if dist ≤ k
        return j, dist
    end

    while (forward && j ≤ min(stop, n)) || (!forward && j ≥ max(stop, 1))
        Eq = Pcom[UInt8(seq[j])+1]
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
function alignback!(H, pat, seq, dist, start, matchstop, forward)
    m = length(pat)
    n = length(seq)

    # initialize the cost column
    for i in 0:m
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
