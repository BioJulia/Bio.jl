# Interval
# ========
#
# Base interval types and utilities.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# Note, just to be clear: this shadows IntervalTrees.Interval
"A genomic interval specifies interval with some associated metadata"
type Interval{T} <: IntervalTrees.AbstractInterval{Int64}
    seqname::StringField
    first::Int64
    last::Int64
    strand::Strand
    metadata::T

    function Interval(seqname, first, last, strand, metadata)
        return new(seqname, first, last, strand, metadata)
    end

    function Interval()
        return new(StringField(), 0, 0, STRAND_NA, T())
    end
end

function Interval{T}(seqname::AbstractString, first::Integer, last::Integer,
                    strand::Union{Strand,Char}, metadata::T)
    return Interval{T}(convert(StringField, seqname), first, last, strand, metadata)
end

function Interval(seqname::AbstractString, first::Integer, last::Integer,
                  strand::Union{Strand,Char}=STRAND_BOTH)
    return Interval{Void}(seqname, first, last, strand, nothing)
end

function Base.copy{T}(interval::Interval{T})
    return Interval{T}(copy(interval.seqname), interval.first, interval.last,
                       interval.strand, copy(interval.metadata))
end

function Bio.seqname(i::Interval)
    return i.seqname
end

function Bio.metadata(i::Interval)
    return i.metadata
end

function strand(i::Interval)
    return i.strand
end

"""
    leftposition(i::Interval)

Return the leftmost position of `i`.
"""
function Bio.leftposition(i::Interval)
    return i.first
end

"""
    rightposition(i::Interval)

Return the rightmost position of `i`.
"""
function Bio.rightposition(i::Interval)
    return i.last
end

IntervalTrees.first(i::Interval) = i.first
IntervalTrees.last(i::Interval) = i.last

function Base.isless{T}(a::Interval{T}, b::Interval{T},
                        seqname_isless::Function=isless)
    if a.seqname != b.seqname
        return seqname_isless(a.seqname, b.seqname)::Bool
    elseif a.first != b.first
        return a.first < b.first
    elseif a.last != b.last
        return a.last < b.last
    elseif a.strand != b.strand
        return a.strand < b.strand
    else
        return false
    end
end

"""
Check if two intervals are well ordered.

Intervals are considered well ordered if a.seqname <= b.seqnamend and
a.first <= b.first.
"""
function isordered{T}(a::Interval{T}, b::Interval{T},
                      seqname_isless::Function=isless)
    if a.seqname != b.seqname
        return seqname_isless(a.seqname, b.seqname)::Bool
    elseif a.first != b.first
        return a.first < b.first
    else
        return true
    end
end

"""
Return true if interval `a` entirely precedes `b`.
"""
function precedes{T}(a::Interval{T}, b::Interval{T},
                     seqname_isless::Function=isless)
    return (a.last < b.first && a.seqname == b.seqname) ||
        seqname_isless(a.seqname, b.seqname)::Bool
end

function Base.:(==){T}(a::Interval{T}, b::Interval{T})
    return a.seqname  == b.seqname &&
           a.first    == b.first &&
           a.last     == b.last &&
           a.strand   == b.strand &&
           a.metadata == b.metadata
end

"Return true if interval `a` overlaps interval `b`, with no consideration to strand"
function Bio.isoverlapping{S, T}(a::Interval{S}, b::Interval{T})
    return a.first <= b.last && b.first <= a.last && a.seqname == b.seqname
end

function Base.show(io::IO, i::Interval)
    if get(io, :compact, false)
        print(io, i.seqname, ":", i.first, "-", i.last, "  ", i.strand, "  ", i.metadata)
    else
        println(io, summary(i), ':')
        println(io, "  sequence name: ", i.seqname)
        println(io, "  leftmost position: ", i.first)
        println(io, "  rightmost position: ", i.last)
        println(io, "  strand: ", i.strand)
          print(io, "  metadata: ", i.metadata)
    end
end

"""
A comparison function used to sort on numbers within text.

This is useful since sequences are often named things like "chr12" or
"read1234". Treating the numbers as numbers and not text gives a more natural
ordering.

This is similar to the '--version-sort' option in GNU coreutils sort.
"""
function alphanum_isless(a::AbstractString, b::AbstractString)
    i = 1
    j = 1

    # match up to the first digit
    k0 = 0 # position of first digit
    while i <= length(a) && j <= length(b)
        if isdigit(a[i]) && isdigit(b[j])
            k0 = i
            break
        else
            if a[i] != b[j]
                return a[i] < b[j]
            end
        end
        i = nextind(a, i)
        j = nextind(b, j)
    end

    # match numbers
    ka1, kb1 = 0, 0
    while i <= length(a) && isdigit(a[i])
        ka1 = i
        i = nextind(a, i)
    end
    while j <= length(b) && isdigit(b[j])
        kb1 = j
        j = nextind(b, j)
    end

    if ka1 == 0 && kb1 != 0
        return true
    elseif ka1 != 0 && kb1 == 0
        return false
    elseif ka1 != 0 && kb1 != 0
        aval = parse(Int, a[k0:ka1])
        bval = parse(Int, b[k0:kb1])
        if aval != bval
            return aval < bval
        end
    end

    # match suffixes
    while i <= length(a) && j <= length(b)
        if a[i] != b[j]
            return a[i] < b[j]
        end
        i = nextind(a, i)
        j = nextind(b, j)
    end

    return j <= length(b)
end

function metadatatype{T}(::Type{T})
    return _metadatatype(eltype(T))
end

function metadatatype(x::Any)
    return metadatatype(typeof(x))
end

function _metadatatype{T}(::Type{Interval{T}})
    return T
end

function _metadatatype{T}(::Type{T})
    return T
end
