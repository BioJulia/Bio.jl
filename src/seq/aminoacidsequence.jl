# Amino acid Sequence
# ===================


# A general purpose amino acid representation.
#
# Amino acid are simple byte arrays using the encoding defined above. Any
# byte outside the range 0x00:0x14 is considered invalid and must result in an
# error.
#
# Like NucleotideSequence, amino acid sequences are immutable by convention.

type AminoAcidSequence
    data::Vector{AminoAcid}
    part::UnitRange{Int} # interval within `data` defining the (sub)sequence
end


# Constructors
# ------------

# Construct a subsequence of another amino acid sequence
function AminoAcidSequence(other::AminoAcidSequence, part::UnitRange)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    if start < other.part.start || stop > other.part.stop
        error("Invalid subsequence range")
    end
    return AminoAcidSequence(other.data, part)
end

# Construct of a subsequence from another amino acid sequence
function AminoAcidSequence(seq::String)
    len = length(seq)
    data = Array(AminoAcid, len)
    for (i, c) in enumerate(seq)
        data[i] = convert(AminoAcid, c)
    end

    return AminoAcidSequence(data, 1:len)
end


# Conversion from/to String
# -------------------------
convert(::Type{AminoAcidSequence}, seq::String) = AminoAcidSequence(seq)
convert(::Type{String}, seq::AminoAcidSequence) = convert(String, [convert(Char, x) for x in seq])


# Basic functions
# ---------------

function show(io::IO, seq::AminoAcidSequence)
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

length(seq::AminoAcidSequence) = length(seq.part)
endof(seq::AminoAcidSequence)  = length(seq)

function getindex(seq::AminoAcidSequence, i::Integer)
    if i > length(seq) || i < 1
        error(BoundsError())
    end
    i += seq.part.start - 1
    return seq.data[i]
end

# Construct a subsequence
getindex(seq::AminoAcidSequence, r::UnitRange) = AminoAcidSequence(seq, r)

# Replace a AminoSequence's data with a copy, copying only what's needed.
function orphan!(seq::AminoAcidSequence, reorphan=true)
    seq.data = seq.data[seq.part]
    seq.part = 1:length(seq.part)
    return seq
end

copy(seq::AminoAcidSequence) = orphan!(AminoAcidSequence(seq.data, seq.part))


# Iterating through amino acid sequence
# -------------------------------------

start(seq::AminoAcidSequence) = seq.part.start

function next(seq::AminoAcidSequence, i)
    aa = seq.data[i]
    return (aa, i + 1)
end

done(seq::AminoAcidSequence, i) = (i > seq.part.stop)


# String decorator
# ----------------

# Enable building sequence literals like: aa"ACDEFMN"
macro aa_str(seq, flags...)
    return AminoAcidSequence(seq)
end
