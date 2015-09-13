
# Base interval types and utilities
# ---------------------------------
bitstype 8 Strand

convert(::Type{Strand}, strand::Uint8) = box(Strand, unbox(Uint8, strand))
convert(::Type{Uint8}, strand::Strand) = box(Uint8, unbox(Strand, strand))


const STRAND_NA   = convert(Strand, 0b000)
const STRAND_POS  = convert(Strand, 0b001)
const STRAND_NEG  = convert(Strand, 0b010)
const STRAND_BOTH = convert(Strand, 0b011)

function show(io::IO, strand::Strand)
    if strand == STRAND_NA
        print(io, "?")
    elseif strand == STRAND_POS
        print(io, "+")
    elseif strand == STRAND_NEG
        print(io, "-")
    elseif strand == STRAND_BOTH
        print(io, ".")
    else
        print(io, "(undefined strand)")
    end
end


isless(a::Strand, b::Strand) = convert(Uint8, a) < convert(Uint8, b)


function convert(::Type{Strand}, strand::Char)
    if strand == '+'
        return STRAND_POS
    elseif strand == '-'
        return STRAND_NEG
    elseif strand == '.'
        return STRAND_BOTH
    elseif strand == '?'
        return STRAND_NA
    else
        error("$(strand) is not a valid strand")
    end
end


# Note, just to be clear: this shadows IntervalTrees.Interval
"A genomic interval specifies interval with some associated metadata"
type Interval{T} <: AbstractInterval{Int64}
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


function Interval{T}(seqname::String, first::Integer, last::Integer,
                    strand::Strand, metadata::T)
    return Interval{T}(convert(StringField, seqname), first, last, strand, metadata)
end


function Interval(seqname::String, first::Integer, last::Integer,
                  strand::Strand=STRAND_BOTH)
    return Interval{Nothing}(seqname, first, last, strand, nothing)
end


function Base.copy{T}(interval::Interval{T})
    return Interval{T}(copy(interval.seqname), interval.first, interval.last,
                       interval.strand, copy(interval.metadata))
end


function first(i::Interval)
    return i.first
end


function last(i::Interval)
    return i.last
end


function isless{T}(a::Interval{T}, b::Interval{T},
                   seqname_isless::Function=alphanum_isless)
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
                      seqname_isless::Function=alphanum_isless)
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
                     seqname_isless::Function=alphanum_isless)
    return (a.last < b.first && a.seqname == b.seqname) ||
        seqname_isless(a.seqname, b.seqname)::Bool
end


function =={T}(a::Interval{T}, b::Interval{T})
    return a.seqname  == b.seqname &&
           a.first    == b.first &&
           a.last     == b.last &&
           a.strand   == b.strand &&
           a.metadata == b.metadata
end


"Return true if interval `a` overlaps interval `b`, with no consideration to strand"
function isoverlapping{S, T}(a::Interval{S}, b::Interval{T})
    return a.first <= b.last && b.first <= a.last && a.seqname == b.seqname
end


function show(io::IO, i::Interval)
    print(io, i.seqname, ":", i.first, "-", i.last, "    ", i.strand, "    ", i.metadata)
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


"""
A type deriving `IntervalStream{T}` must be iterable and produce
Interval{T} objects in sorted order.
"""
abstract IntervalStream{T}


typealias IntervalStreamOrArray{T} Union(Vector{Interval{T}}, IntervalStream{T},
                                         AbstractParser)


function metadatatype{T}(::IntervalStream{T})
    return T
end
