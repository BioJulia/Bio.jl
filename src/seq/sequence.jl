# Sequence
# ========
#
# Abstract biological sequence type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
Abstract biological sequence type.

Any subtype `S <: Sequence` should implement the following methods:

* `Base.length(seq::S)`: return the length of `seq`.
* `Base.eltype(::Type{S})`: return the element type of `S`.
* `inbounds_getindex(seq::S, i::Integer)`: return the element at `i` of `seq`
    without checking bounds
"""
abstract Sequence

# This is useful for obscure reasons. We use SeqRecord{Sequence} for reading
# sequence in an undetermined alphabet, but a consequence that we need to be
# able to construct a `Sequence`.
function Sequence()
    return DNASequence()
end

Base.size(seq::Sequence) = (length(seq),)
Base.endof(seq::Sequence) = length(seq)
Base.isempty(seq::Sequence) = length(seq) == 0
Base.eachindex(seq::Sequence) = 1:endof(seq)

@inline function Base.checkbounds(seq::Sequence, i::Integer)
    if 1 ≤ i ≤ endof(seq)
        return true
    end
    throw(BoundsError(seq, i))
end

if VERSION < v"0.5-"
    macro boundscheck(ex)
        # no op
        return ex
    end
end

function Base.getindex(seq::Sequence, i::Integer)
    @boundscheck checkbounds(seq, i)
    return inbounds_getindex(seq, i)
end

Base.start(seq::Sequence) = 1
Base.done(seq::Sequence, i) = i > endof(seq)
Base.next(seq::Sequence, i) = inbounds_getindex(seq, i), i + 1


# Comparison
# ----------

@compat function Base.:(==)(seq1::Sequence, seq2::Sequence)
    return eltype(seq1)    == eltype(seq2) &&
           length(seq1)    == length(seq2) &&
           cmp(seq1, seq2) == 0
end

Base.isless(seq1::Sequence, seq2::Sequence) = cmp(seq1, seq2) < 0

function Base.cmp(seq1::Sequence, seq2::Sequence)
    m = length(seq1)
    n = length(seq2)
    for i in 1:min(m, n)
        c = cmp(inbounds_getindex(seq1, i),
                inbounds_getindex(seq2, i))
        if c != 0
            return c
        end
    end
    return cmp(m, n)
end


# Finders
# -------

function Base.findnext(seq::Sequence, val, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:endof(seq)
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return 0
end

function Base.findprev(seq::Sequence, val, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:-1:1
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return 0
end

Base.findfirst(seq::Sequence, val) = findnext(seq, val, 1)
Base.findlast(seq::Sequence, val) = findprev(seq, val, endof(seq))


"""
    gc_content(seq::Sequence)

Calculate GC content of `seq`.
"""
function gc_content(seq::Sequence)
    if !(eltype(seq) <: Nucleotide)
        error("elements must be nucleotides")
    end

    gc = 0
    for x in seq
        if isGC(x)
            gc += 1
        end
    end

    if isempty(seq)
        return 0.0
    else
        return gc / length(seq)
    end
end

# Printers
# --------

function Base.print(io::IO, seq::Sequence)
    width = 70
    col = 1
    for x in seq
        if col > width
            write(io, '\n')
            col = 1
        end
        print(io, x)
        col += 1
    end
end

function Base.show(io::IO, seq::Sequence)
    println(io, summary(seq), ':')
    # don't show more than this many characters
    # to avoid filling the screen with junk
    width = displaysize()[2]
    if length(seq) > width
        half = div(width, 2)
        for i in 1:half-1
            print(io, seq[i])
        end
        print(io, '…')
        for i in endof(seq)-half+2:endof(seq)
            print(io, seq[i])
        end
    else
        for x in seq
            print(io, Char(x))
        end
    end
end
