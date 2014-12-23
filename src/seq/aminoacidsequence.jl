# Aminoacid Sequence
# ==================


# A general purpose amino acid representation.
#
# Amino acid are simple byte arrays using the encoding defined above. Any
# byte outside the range 0x00:0x14 is considered invalid and must result in an
# error.
#
# Like NucleotideSequence, amino acid sequences are immutable by convention.

type AminoacidSequence
    data::Vector{AminoAcid}
    part::UnitRange{Int} # interval within `data` defining the (sub)sequence
end


# Constructors
# ------------

# Construct a subsequence of another aminoacid sequence
function AminoacidSequence(other::AminoacidSequence, part::UnitRange)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    if start < other.part.start || stop > other.part.stop
        error("Invalid subsequence range")
    end
    return AminoacidSequence(other.data, part)
end

# Construct of a subsequence from another amino acid sequence
function AminoacidSequence(seq::String)
    len = length(seq)
    data = Array(AminoAcid, len)
    for (i, c) in enumerate(seq)
        data[i] = convert(AminoAcid, c)
    end

    return AminoacidSequence(data, 1:len)
end


# Conversion from/to String
# -------------------------
convert(::Type{AminoacidSequence}, seq::String) = AminoacidSequence(seq)
convert(::Type{String}, seq::AminoacidSequence) = convert(String, [convert(Char, x) for x in seq])


# Basic functions
# ---------------

function show(io::IO, seq::AminoacidSequence)
    len = length(seq)
    write(io, "$(string(len))aa Sequence:\n ")

    const maxcount = 50
    if len > maxcount
        for aa in seq[1:div(maxcount, 2) - 1]
            write(io, convert(Char, aa))
        end
        write(io, "…")
        for aa in seq[(end - (div(maxcount, 2) - 1)):end]
            write(io, convert(Char, aa))
        end
    else
        for aa in seq
            write(io, convert(Char, aa))
        end
    end
end

length(seq::AminoacidSequence) = length(seq.part)
endof(seq::AminoacidSequence)  = length(seq)

function getindex(seq::AminoacidSequence, i::Integer)
    if i > length(seq) || i < 1
        error(BoundsError())
    end
    i += seq.part.start - 1
    return seq.data[i]
end

# Construct a subsequence
getindex(seq::AminoacidSequence, r::UnitRange) = AminoacidSequence(seq, r)

# Replace a AminoSequence's data with a copy, copying only what's needed.
function orphan!(seq::AminoacidSequence, reorphan=true)
    seq.data = seq.data[seq.part]
    seq.part = 1:length(seq.part)
    return seq
end

copy(seq::AminoacidSequence) = orphan!(AminoacidSequence(seq.data, seq.part))


# Iterating through aminoacid sequence
# ------------------------------------

start(seq::AminoacidSequence) = seq.part.start

function next(seq::AminoacidSequence, i)
    aa = seq.data[i]
    return (aa, i + 1)
end

done(seq::AminoacidSequence, i) = (i > seq.part.stop)


# String decorator
# ----------------

# Enable building sequence literals like: aa"ACDEFMN"
macro aa_str(seq, flags...)
    return AminoacidSequence(seq)
end
