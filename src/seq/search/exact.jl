# Exact Search
# ============

# Reference Implementation
# ------------------------
#
# `searchindex` and `rsearchindex` should work exactly the same way as the
# following reference implementation:
#=
# Return the index of the first occurrence of `pat` in `seq[start:stop]`.
function searchindex(seq, pat, start=1, stop=endof(seq))
    m = length(pat)
    n = length(seq)
    for s in max(start-1, 0):min(stop, n)-m
        if occurs_with_shift(pat, seq, s)
            return s+1  # found
        end
    end
    return 0  # not found
end

# Return the index of the last occurrence of `pat` in `seq[stop:start]`.
function rsearchindex(seq, pat, start=endof(seq), stop=1)
    n = length(seq)
    m = length(pat)
    for s in min(start, n)-m:-1:max(stop-1, 0)
        if occurs_with_shift(pat, seq, s)
            return s+1  # found
        end
    end
    return 0  # not found
end

function occurs_with_shift(pat, seq, shift)
    for i in 1:endof(pat)
        if !iscompatible(pat[i], seq[shift+i])
            return false
        end
    end
    return true
end
=#


# Actual Implementation
# ---------------------

"""
Query type for exact sequence search.
"""
immutable ExactSearchQuery{S<:Sequence}
    seq::S         # query sequence
    cbits::UInt32  # compatibility bits
    fshift::Int    # shift length for forward search
    bshift::Int    # shift length for backward search
end

function ExactSearchQuery(query::Sequence)
    cbits, fshift, bshift = preprocess(query)
    return ExactSearchQuery(query, cbits, fshift, bshift)
end

function preprocess(query)
    if length(query) == 0
        return UInt32(0), 0, 0
    end
    m = length(query)
    first = query[1]
    last = query[end]
    cbits::UInt32 = 0
    fshift = bshift = m
    for i in 1:endof(query)
        x = query[i]
        cbits |= compatbits(x)
        if iscompatible(x, last) && i < m
            fshift = m - i
        end
    end
    for i in endof(query):-1:1
        x = query[i]
        if iscompatible(x, first) && i > 1
            bshift = i - 1
        end
    end
    return cbits, fshift, bshift
end

function checkeltype(seq1, seq2)
    if eltype(seq1) != eltype(seq2)
        throw(ArgumentError("the element type of two sequences must match"))
    end
end


# Forward
# -------


# NOTE: the return value is a range while the counterpart for strings is an index
"""
    search(seq::Sequence, pat[, start=1[, stop=endof(seq)]])

Return the range of the first occurrence of `pat` in `seq[start:stop]`; symbol
comparison is done using `Bio.Seq.iscompatible`.
"""
function Base.search(seq::Sequence, val,
                     start::Integer=1, stop::Integer=endof(seq))
    s = searchindex(seq, val, start, stop)
    if s == 0
        return 0:-1
    else
        return s:s
    end
end

function Base.search(seq::Sequence, pat::Sequence,
                     start::Integer=1, stop::Integer=endof(seq))
    return search(seq, ExactSearchQuery(pat), start, stop)
end

function Base.search(seq::Sequence, query::ExactSearchQuery,
                     start::Integer=1, stop::Integer=endof(seq))
    s = searchindex(seq, query, start, stop)
    if s == 0
        return 0:-1
    end
    return s:s+length(query.seq)-1
end

"""
    searchindex(seq::Sequence, pat[, start=1[, stop=endof(seq)]])

Return the index of the first occurrence of `pat` in `seq[start:stop]`; symbol
comparison is done using `Bio.Seq.iscompatible`.
"""
function Base.searchindex(seq::Sequence, val,
                          start::Integer=1, stop::Integer=endof(seq))
    x = convert(eltype(seq), val)
    for i in max(start, 1):min(stop, endof(seq))
        if iscompatible(seq[i], x)
            return i  # found
        end
    end
    return 0  # not found
end

function Base.searchindex(seq::Sequence, pat::Sequence,
                          start::Integer=1, stop::Integer=endof(seq))
    return searchindex(seq, ExactSearchQuery(pat), start, stop)
end

function Base.searchindex(seq::Sequence, query::ExactSearchQuery,
                          start::Integer=1, stop::Integer=endof(seq))
    checkeltype(seq, query.seq)
    return quicksearch(query, seq, start, stop)
end

# This algorithm is borrowed from base/strings/search.jl, which looks similar to
# Sunday's Quick Search algorithm.
function quicksearch(query, seq, start, stop)
    pat = query.seq
    m = length(pat)
    n = length(seq)
    stop′ = min(stop, n) - m
    s::Int = max(start - 1, 0)

    if m == 0  # empty query
        if s ≤ stop′
            return s + 1  # found
        else
            return 0  # not found
        end
    end

    while s ≤ stop′
        if iscompatible(pat[m], seq[s+m])
            i = m - 1
            while i > 0
                if !iscompatible(pat[i], seq[s+i])
                    break
                end
                i -= 1
            end
            if i == 0
                return s + 1  # found
            elseif s < stop′ && (query.cbits & compatbits(seq[s+m+1]) == 0)
                s += m + 1
            elseif isambiguous(seq[s+m])
                s += 1
            else
                s += query.fshift
            end
        elseif s < stop′ && (query.cbits & compatbits(seq[s+m+1]) == 0)
            s += m + 1
        else
            s += 1
        end
    end

    return 0  # not found
end


# Backward
# --------

# NOTE: the return value is a range while the counterpart for strings is an index
"""
    rsearch(seq::Sequence, pat[, start=endof(seq)[, stop=1]])

Return the range of the last occurrence of `pat` in `seq[start:stop]`; symbol
comparison is done using `Bio.Seq.iscompatible`.
"""
function Base.rsearch(seq::Sequence, val,
                      start::Integer=endof(seq), stop::Integer=1)
    s = rsearchindex(seq, val, start, stop)
    if s == 0
        return 0:-1
    else
        return s:s
    end
end

function Base.rsearch(seq::Sequence, pat::Sequence,
                      start::Integer=endof(seq), stop::Integer=1)
    return rsearch(seq, ExactSearchQuery(pat), start, stop)
end

function Base.rsearch(seq::Sequence, query::ExactSearchQuery,
                      start::Integer=endof(seq), stop::Integer=1)
    s = rsearchindex(seq, query, start, stop)
    if s == 0
        return 0:-1
    end
    return s:s+length(query.seq)-1
end

"""
    rsearchindex(seq::Sequence, pat[, start=1[, stop=endof(seq)]])

Return the index of the last occurrence of `pat` in `seq[start:stop]`; symbol
comparison is done using `Bio.Seq.iscompatible`.
"""
function Base.rsearchindex(seq::Sequence, val,
                           start::Integer=endof(seq), stop::Integer=1)
    x = convert(eltype(seq), val)
    for i in min(start, endof(seq)):-1:max(stop, 1)
        if iscompatible(seq[i], x)
            return i  # found
        end
    end
    return 0  # not found
end

function Base.rsearchindex(seq::Sequence, pat::Sequence,
                           start::Integer=endof(seq), stop::Integer=1)
    return rsearchindex(seq, ExactSearchQuery(pat), start, stop)
end

function Base.rsearchindex(seq::Sequence, query::ExactSearchQuery,
                           start::Integer=endof(seq), stop::Integer=1)
    checkeltype(seq, query.seq)
    return quickrsearch(seq, query, start, stop)
end

function quickrsearch(seq, query, start, stop)
    pat = query.seq
    m = length(pat)
    n = length(seq)
    stop′ = max(stop - 1, 0)
    s::Int = min(start, n) - m

    if m == 0  # empty query
        if s ≥ stop′
            return s + 1  # found
        else
            return 0  # not found
        end
    end

    while s ≥ stop′
        if iscompatible(pat[1], seq[s+1])
            i = 2
            while i < m + 1
                if !iscompatible(pat[i], seq[s+i])
                    break
                end
                i += 1
            end
            if i == m + 1
                return s + 1  # found
            elseif s > stop′ && (query.cbits & compatbits(seq[s]) == 0)
                s -= m + 1
            elseif isambiguous(seq[s+1])
                s -= 1
            else
                s -= query.bshift
            end
        elseif s > stop′ && (query.cbits & compatbits(seq[s]) == 0)
            s -= m + 1
        else
            s -= 1
        end
    end

    return 0  # not found
end
