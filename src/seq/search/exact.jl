# Exact Search
# ============

function checkeltype(seq1, seq2)
    if eltype(seq1) != eltype(seq2)
        throw(ArgumentError("the element type of two sequences must match"))
    end
    return true
end


# Forward
# -------

function Base.search(seq::Sequence, val, start::Integer=1)
    v = convert(eltype(seq), val)
    for i in Int(max(start,1)):endof(seq)
        if iscompatible(seq[i], v)
            return i
        end
    end
    return 0
end

function Base.search(seq::Sequence, pat::Sequence, start::Integer=1)
    s = searchindex(seq, pat, start)
    if s == 0
        return 0:-1
    end
    return s:s+length(pat)-1
end

function Base.searchindex(seq::Sequence, pat::Sequence, start::Integer=1)
    checkeltype(seq, pat)
    return quicksearch(seq, pat, start)
end

# This algorithm is borrowed from base/strings/search.jl, which looks similar to
# Sunday's Quick Search algorithm.
function quicksearch(seq, pat, start)
    m = length(pat)
    n = length(seq)

    if m == 0
        return Int(start)
    end

    # preprocess
    last = pat[m]
    cbits::UInt32 = 0
    shift = m
    for i in 1:m
        x = pat[i]
        cbits |= compatbits(x)
        if iscompatible(x, last) && i < m
            shift = m - i
        end
    end

    j::Int = max(start, 1)
    @inbounds while j ≤ n - m + 1
        if iscompatible(seq[j+m-1], last)
            i = 1
            while i ≤ m - 1
                if !iscompatible(seq[j+i-1], pat[i])
                    break
                end
                i += 1
            end

            if i == m
                # found
                return j
            elseif j ≤ n - m && (cbits & compatbits(seq[j+m]) == 0)
                j += m + 1
            elseif isambiguous(seq[j+m-1])
                j += 1
            else
                j += shift
            end
        elseif j ≤ n - m && (cbits & compatbits(seq[j+m]) == 0)
            j += m + 1
        else
            j += 1
        end
    end

    # not found
    return 0
end


# Backward
# --------

function Base.rsearch(seq::Sequence, val, start::Integer=endof(seq))
    v = convert(eltype(seq), val)
    for i in Int(min(start,endof(seq))):-1:1
        if iscompatible(seq[i], v)
            return i
        end
    end
    return 0
end

function Base.rsearch(seq::Sequence, pat::Sequence, start::Integer=endof(seq))
    s = rsearchindex(seq, pat, start)
    if s == 0
        return 0:-1
    end
    return s:s+length(pat)-1
end

function Base.rsearchindex(seq::Sequence, pat::Sequence, start::Integer=endof(seq))
    checkeltype(seq, pat)
    return quickrsearch(seq, pat, start)
end

function quickrsearch(seq, pat, start)
    m = length(pat)
    n = length(seq)

    if m == 0
        return Int(start)
    end

    # preprocess
    first = pat[1]
    cbits::UInt32 = 0
    shift = m
    for i in m:-1:1
        x = pat[i]
        cbits |= compatbits(x)
        if iscompatible(x, first) && i > 1
            shift = i - 1
        end
    end

    j::Int = min(start, endof(seq))
    @inbounds while j ≥ m
        if iscompatible(seq[j-m+1], first)
            i = m
            while i > 1
                if !iscompatible(seq[j+i-m], pat[i])
                    break
                end
                i -= 1
            end

            if i == 1
                # found
                return j - m + 1
            elseif j > m && (cbits & compatbits(seq[j-m]) == 0)
                j -= m + 1
            elseif isambiguous(seq[j-m+1])
                j -= 1
            else
                j -= shift
            end
        elseif j > m && (cbits & compatbits(seq[j-m]) == 0)
            j -= m + 1
        else
            j -= 1
        end
    end

    # not found
    return 0
end
