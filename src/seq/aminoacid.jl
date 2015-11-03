# Aminoacids
# ==========

"Type representing AminoAcids"
bitstype 8 AminoAcid


# Conversion from/to integers
# ---------------------------

convert(::Type{AminoAcid}, aa::UInt8) = box(AminoAcid, unbox(UInt8, aa))
convert(::Type{UInt8}, aa::AminoAcid) = box(UInt8, unbox(AminoAcid, aa))
convert{T <: Unsigned}(::Type{T}, aa::AminoAcid) = box(T, Base.zext_int(T, unbox(AminoAcid, aa)))
convert{T <: Unsigned}(::Type{AminoAcid}, aa::T) = convert(AminoAcid, convert(UInt8, aa))


# Amino acid encoding definition
# ------------------------------


"Alanine"
const AA_A = convert(AminoAcid, 0x00)

"Arginine"
const AA_R = convert(AminoAcid, 0x01)

"Asparagine"
const AA_N = convert(AminoAcid, 0x02)

"Aspartic Acid"
const AA_D = convert(AminoAcid, 0x03)

"Cysteine"
const AA_C = convert(AminoAcid, 0x04)

"Glutamine"
const AA_Q = convert(AminoAcid, 0x05)

"Glutamic Acid"
const AA_E = convert(AminoAcid, 0x06)

"Glycine"
const AA_G = convert(AminoAcid, 0x07)

"Histidine"
const AA_H = convert(AminoAcid, 0x08)

"Isoleucine"
const AA_I = convert(AminoAcid, 0x09)

"Leucine"
const AA_L = convert(AminoAcid, 0x0a)

"Lysine"
const AA_K = convert(AminoAcid, 0x0b)

"Methionine"
const AA_M = convert(AminoAcid, 0x0c)

"Phenylalanine"
const AA_F = convert(AminoAcid, 0x0d)

"Proline"
const AA_P = convert(AminoAcid, 0x0e)

"Serine"
const AA_S = convert(AminoAcid, 0x0f)

"Threonine"
const AA_T = convert(AminoAcid, 0x10)

"Tryptophan"
const AA_W = convert(AminoAcid, 0x11)

"Tyrosine"
const AA_Y = convert(AminoAcid, 0x12)

"Valine"
const AA_V = convert(AminoAcid, 0x13)

"Pyrrolysine"
const AA_O = convert(AminoAcid, 0x14)

"Selenocysteine"
const AA_U = convert(AminoAcid, 0x15)

"Aspartic Acid or Asparagine"
const AA_B = convert(AminoAcid, 0x16)

"Glutamine or Glutamic Acid"
const AA_Z = convert(AminoAcid, 0x17)

"Unspecified or Unknown Amino Acid"
const AA_X = convert(AminoAcid, 0x18)

"Invalid Amino Acid"
const AA_INVALID = convert(AminoAcid, 0x19) # Used during conversion from strings


function isvalid(aa::AminoAcid)
    return convert(UInt8, aa) ≤ convert(UInt8, AA_X)
end


# Conversion from/to Char
# -----------------------

# lookup table for characters in 'A':'z'
const char_to_aa = [
    AA_A,       AA_B,       AA_C,       AA_D,       AA_E,       AA_F,
    AA_G,       AA_H,       AA_I,       AA_INVALID, AA_K,       AA_L,
    AA_M,       AA_N,       AA_O,       AA_P,       AA_Q,       AA_R,
    AA_S,       AA_T,       AA_U,       AA_V,       AA_W,       AA_X,
    AA_Y,       AA_Z,       AA_INVALID, AA_INVALID, AA_INVALID, AA_INVALID,
    AA_INVALID, AA_INVALID, AA_A,       AA_B,       AA_C,       AA_D,
    AA_E,       AA_F,       AA_G,       AA_H,       AA_I,       AA_INVALID,
    AA_K,       AA_L,       AA_M,       AA_N,       AA_O,       AA_P,
    AA_Q,       AA_R,       AA_S,       AA_T,       AA_U,       AA_V,
    AA_W,       AA_X,       AA_Y,       AA_Z]

function convert(::Type{AminoAcid}, c::Char)
    @inbounds aa = 'A' <= c <= 'y' ? char_to_aa[c - 'A' + 1] : AA_INVALID
    @assert aa != AA_INVALID error("$(c) is not a valid amino acid")
    return aa
end

const aa_to_char = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'O', 'U',
    'B', 'Z', 'X' ]
convert(::Type{Char}, aa::AminoAcid) = aa_to_char[convert(UInt8, aa) + 1]


# Basic functions
# ---------------

function show(io::IO, aa::AminoAcid)
    if aa == AA_INVALID
        write(io, "Invalid Amino Acid")
    else
        write(io, convert(Char, aa))
    end
end

# lookup table of 20 standard amino acids
const threeletter_to_aa = Dict(
    "ALA" => AA_A, "ARG" => AA_R, "ASN" => AA_N, "ASP" => AA_D, "CYS" => AA_C,
    "GLN" => AA_Q, "GLU" => AA_E, "GLY" => AA_G, "HIS" => AA_H, "ILE" => AA_I,
    "LEU" => AA_L, "LYS" => AA_K, "MET" => AA_M, "PHE" => AA_F, "PRO" => AA_P,
    "SER" => AA_S, "THR" => AA_T, "TRP" => AA_W, "TYR" => AA_Y, "VAL" => AA_V, "PYL" => AA_O, "SEC" => AA_U,
    "ASX" => AA_B, "GLX" => AA_Z,
)

function parse(::Type{AminoAcid}, s::AbstractString)
    s′ = strip(s)
    if length(s′) == 1
        return convert(AminoAcid, s′[1])
    end
    try
        return threeletter_to_aa[uppercase(s′)]
    catch ex
        if isa(ex, KeyError)
            error("invalid amino acid string: \"$s\" ")
        end
        rethrow()
    end
end

# TODO: tryparse

# Amino Acids Sequences
# =====================


# A general purpose amino acid representation.
#
# Amino acid are simple byte arrays using the encoding defined above. Any
# byte outside the range 0x00:0x14 is considered invalid and must result in an
# error.
#
# Like NucleotideSequence, amino acid sequences are immutable by convention.

"Type representing AminoAcid Sequences"
type AminoAcidSequence <: Sequence
    data::Vector{AminoAcid}
    part::UnitRange{Int} # interval within `data` defining the (sub)sequence
    mutable::Bool # true if the sequence can be safely mutated

    # true if this was constructed as a subsequence of another sequence or if
    # subsequences were constructed from this sequence. When this is true, we
    # need to copy the data to convert from immutable to mutable
    hasrelatives::Bool
end


# Constructors
# ------------

"Construct a subsequence of another amino acid sequence"
function AminoAcidSequence(other::AminoAcidSequence, part::UnitRange;
                           mutable::Bool=false)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    if start < other.part.start || stop > other.part.stop
        error("Invalid subsequence range")
    end
    seq = AminoAcidSequence(other.data, part, mutable, true)

    if other.mutable || mutable
        orphan!(seq, true)
    end

    other.hasrelatives = !other.mutable
    seq.hasrelatives = !mutable

    return seq
end


"Construct of a subsequence from another amino acid sequence"
function AminoAcidSequence(seq::Union{Vector{UInt8}, AbstractString},
                           startpos::Int, endpos::Int, unsafe::Bool=false;
                           mutable::Bool=false)

    len = endpos - startpos + 1
    data = Array(AminoAcid, len)
    for (i, j) in enumerate(startpos:endpos)
        data[i] = convert(AminoAcid, convert(Char, seq[j]))
    end

    return AminoAcidSequence(data, 1:len, mutable, false)
end

function AminoAcidSequence()
    return AminoAcidSequence(AminoAcid[], 1:0, true, false)
end

function AminoAcidSequence(seq::AbstractVector{AminoAcid},
                           startpos::Int, endpos::Int, unsafe::Bool=false;
                           mutable::Bool=false)
    len = endpos - startpos + 1
    data = Vector{AminoAcid}(len)
    for (i, j) in enumerate(startpos:endpos)
        aa = seq[j]
        if !unsafe && !isvalid(aa)
            error("the sequence includes an invalid amino acid at $j")
        end
        data[i] = aa
    end
    return AminoAcidSequence(data, 1:len, mutable, false)
end

"Construct an amino acid sequence by concatenating other sequences"
function AminoAcidSequence(chunks::AminoAcidSequence...)
    seqlen = 0
    for chunk in chunks
        seqlen += length(chunk)
    end

    data = Array(AminoAcid, seqlen)
    pos = 1
    for chunk in chunks
        copy!(data, pos, chunk.data, chunk.part.start, length(chunk))
        pos += length(chunk)
    end

    return AminoAcidSequence(data, 1:seqlen, false, false)
end


(*)(chunk1::AminoAcidSequence, chunks::AminoAcidSequence...) = AminoAcidSequence(chunk1, chunks...)


"Construct an amino acid sequence by repeating another sequence"
function repeat(chunk::AminoAcidSequence, n::Integer)
    seqlen = n * length(chunk)

    data = Array(AminoAcid, seqlen)
    pos = 1
    for i in 1:n
        copy!(data, pos, chunk.data, chunk.part.start, length(chunk))
        pos += length(chunk)
    end

    return AminoAcidSequence(data, 1:seqlen, false, false)
end


(^)(chunk::AminoAcidSequence, n::Integer) = repeat(chunk, n::Integer)


# Conversion
# ----------

# Conversion from/to a byte sequence
convert(::Type{AminoAcidSequence}, seq::AbstractVector{AminoAcid}) = AminoAcidSequence(seq, 1, endof(seq))
convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = [convert(AminoAcid, x) for x in seq]

# Conversion from/to String
convert(::Type{AminoAcidSequence}, seq::AbstractString) = AminoAcidSequence(seq, 1, endof(seq))
convert(::Type{AbstractString}, seq::AminoAcidSequence) = convert(AbstractString, [convert(Char, x) for x in seq])


# Basic functions
# ---------------

function ==(a::AminoAcidSequence, b::AminoAcidSequence)
    if (len = length(a)) != length(b)
        return false
    end
    for i in 1:len
        if a[i] != b[i]
            return false
        end
    end
    return true
end


function show(io::IO, seq::AminoAcidSequence)
    len = length(seq)
    write(io, "$(string(len))aa ",
          seq.mutable ? "Mutable " : "",
          "Sequence:\n ")

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
function orphan!(seq::AminoAcidSequence, reorphan=false)
    if !reorphan && seq.part.start == 1
        return seq
    end

    seq.data = seq.data[seq.part]
    seq.part = 1:length(seq.part)
    seq.hasrelatives = false
    return seq
end


copy(seq::AminoAcidSequence) = orphan!(AminoAcidSequence(seq.data, seq.part, seq.mutable, false), true)


function setindex!(seq::AminoAcidSequence, nt::AminoAcid, i::Integer)
    if !seq.mutable
        error("Cannot mutate an immutable sequence. Call `mutable!(seq)` first.")
    end
    seq.data[i] = nt
end


function setindex!(seq::AminoAcidSequence, nt::Char, i::Integer)
    setindex!(seq, convert(AminoAcid, nt), i)
end


"""
Reset the contents of a mutable sequence from a string.
"""
function copy!(seq::AminoAcidSequence, strdata::Vector{UInt8},
               startpos::Integer, stoppos::Integer)
    if !seq.mutable
        error("Cannot copy! to immutable sequnce. Call `mutable!(seq)` first.")
    end

    n = stoppos - startpos + 1
    if length(seq.data) < n
        resize!(seq.data, n)
    end

    for i in 1:n
        seq.data[i] = convert(AminoAcid, Char(strdata[startpos + i - 1]))
    end
    seq.part = 1:n
end


# Mutability/Immutability
# -----------------------

ismutable(seq::AminoAcidSequence) = return seq.mutable


function mutable!(seq::AminoAcidSequence)
    if !seq.mutable
        if seq.hasrelatives
            orphan!(seq, true)
        end
        seq.mutable = true
    end
    return seq
end


function immutable!(seq::AminoAcidSequence)
    seq.mutable = false
    return seq
end


# Iterating through amino acid sequence
# -------------------------------------

start(seq::AminoAcidSequence) = seq.part.start
next(seq::AminoAcidSequence, i) = (seq.data[i], i+1)


function next(seq::AminoAcidSequence, i)
    aa = seq.data[i]
    return (aa, i + 1)
end


done(seq::AminoAcidSequence, i) = (i > seq.part.stop)


# String decorator
# ----------------

# Enable building sequence literals like: aa"ACDEFMN"
macro aa_str(seq, flags...)
    return AminoAcidSequence(remove_newlines(seq))
end
