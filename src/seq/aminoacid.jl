# Aminoacids
# ==========


bitstype 8 AminoAcid


# Conversion from/to integers
# ---------------------------

convert(::Type{AminoAcid}, aa::Uint8) = box(AminoAcid, unbox(Uint8, aa))
convert(::Type{Uint8}, aa::AminoAcid) = box(Uint8, unbox(AminoAcid, aa))
convert{T <: Unsigned}(::Type{T}, aa::AminoAcid) = box(T, Base.zext_int(T, unbox(AminoAcid, aa)))
convert{T <: Unsigned}(::Type{AminoAcid}, aa::T) = convert(AminoAcid, convert(Uint8, aa))


# Amino acid encoding definition
# ------------------------------


@doc doc"Alanine" ->
const AA_A = convert(AminoAcid, 0x00)

@doc doc"Arginine" ->
const AA_R = convert(AminoAcid, 0x01)

@doc doc"Asparagine" ->
const AA_N = convert(AminoAcid, 0x02)

@doc doc"Aspartic Acid" ->
const AA_D = convert(AminoAcid, 0x03)

@doc doc"Cysteine" ->
const AA_C = convert(AminoAcid, 0x04)

@doc doc"Glutamine" ->
const AA_Q = convert(AminoAcid, 0x05)

@doc doc"Glutamic Acid" ->
const AA_E = convert(AminoAcid, 0x06)

@doc doc"Glycine" ->
const AA_G = convert(AminoAcid, 0x07)

@doc doc"Histidine" ->
const AA_H = convert(AminoAcid, 0x08)

@doc doc"Isoleucine" ->
const AA_I = convert(AminoAcid, 0x09)

@doc doc"Leucine" ->
const AA_L = convert(AminoAcid, 0x0a)

@doc doc"Lysine" ->
const AA_K = convert(AminoAcid, 0x0b)

@doc doc"Methionine" ->
const AA_M = convert(AminoAcid, 0x0c)

@doc doc"Phenylalanine" ->
const AA_F = convert(AminoAcid, 0x0d)

@doc doc"Proline" ->
const AA_P = convert(AminoAcid, 0x0e)

@doc doc"Serine" ->
const AA_S = convert(AminoAcid, 0x0f)

@doc doc"Threonine" ->
const AA_T = convert(AminoAcid, 0x10)

@doc doc"Tryptophan" ->
const AA_W = convert(AminoAcid, 0x11)

@doc doc"Tyrosine" ->
const AA_Y = convert(AminoAcid, 0x12)

@doc doc"Valine" ->
const AA_V = convert(AminoAcid, 0x13)

@doc doc"Unspecified or Unknown Amino Acid" ->
const AA_X = convert(AminoAcid, 0x14)

@doc doc"Invalid Amino Acid" ->
const AA_INVALID = convert(AminoAcid, 0x15) # Used during conversion from strings


# Conversion from/to Char
# -----------------------

# lookup table for characters in 'A':'y'
const char_to_aa = [
    AA_A,       AA_INVALID, AA_C,       AA_D,       AA_E,       AA_F,
    AA_G,       AA_H,       AA_I,       AA_INVALID, AA_K,       AA_L,
    AA_M,       AA_N,       AA_INVALID, AA_P,       AA_Q,       AA_R,
    AA_S,       AA_T,       AA_INVALID, AA_V,       AA_W,       AA_X,
    AA_Y,       AA_INVALID, AA_INVALID, AA_INVALID, AA_INVALID, AA_INVALID,
    AA_INVALID, AA_INVALID, AA_A,       AA_INVALID, AA_C,       AA_D,
    AA_E,       AA_F,       AA_G,       AA_H,       AA_I,       AA_INVALID,
    AA_K,       AA_L,       AA_M,       AA_N,       AA_INVALID, AA_P,
    AA_Q,       AA_R,       AA_S,       AA_T,       AA_INVALID, AA_V,
    AA_W,       AA_X,       AA_Y ]

function convert(::Type{AminoAcid}, c::Char)
    @inbounds aa = 'A' <= c <= 'y' ? char_to_aa[c - 'A' + 1] : AA_INVALID
    @assert aa != AA_INVALID error("$(c) is not a valid amino acid")
    return aa
end

const aa_to_char = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X' ]
convert(::Type{Char}, aa::AminoAcid) = aa_to_char[convert(Uint8, aa) + 1]


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
const threeletter_to_aa = @compat Dict(
    "ALA" => AA_A, "ARG" => AA_R, "ASN" => AA_N, "ASP" => AA_D, "CYS" => AA_C,
    "GLN" => AA_Q, "GLU" => AA_E, "GLY" => AA_G, "HIS" => AA_H, "ILE" => AA_I,
    "LEU" => AA_L, "LYS" => AA_K, "MET" => AA_M, "PHE" => AA_F, "PRO" => AA_P,
    "SER" => AA_S, "THR" => AA_T, "TRP" => AA_W, "TYR" => AA_Y, "VAL" => AA_V,
)

function parse(::Type{AminoAcid}, s::String)
    s′ = strip(s)
    if length(s′) == 1
        return convert(AminoAcid, s′[1])
    end
    try
        return threeletter_to_aa[uppercase(s′)]
    catch ex
        if isa(ex, KeyError)
            error("invalid amino acid string: \"$s\"")
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
function AminoAcidSequence(seq::Union(Vector{Uint8}, String),
                           startpos::Int, endpos::Int, unsafe::Bool=false)
    len = endpos - startpos + 1
    data = Array(AminoAcid, len)
    for (i, j) in enumerate(startpos:endpos)
        data[i] = convert(AminoAcid, convert(Char, seq[j]))
    end

    return AminoAcidSequence(data, 1:len)
end

@doc """
Construct an amino acid sequence by concatenating other sequences.
""" ->
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

    return AminoAcidSequence(data, 1:seqlen)
end


(*)(chunk1::AminoAcidSequence, chunks::AminoAcidSequence...) = AminoAcidSequence(chunk1, chunks...)


@doc """
Construct an amino acid sequence by repeating another sequence.
""" ->
function repeat(chunk::AminoAcidSequence, n::Integer)
    seqlen = n * length(chunk)

    data = Array(AminoAcid, seqlen)
    pos = 1
    for i in 1:n
        copy!(data, pos, chunk.data, chunk.part.start, length(chunk))
        pos += length(chunk)
    end

    return AminoAcidSequence(data, 1:seqlen)
end


(^)(chunk::AminoAcidSequence, n::Integer) = repeat(chunk, n::Integer)


function AminoAcidSequence(seq::Union(Vector{Uint8}, String))
    return AminoAcidSequence(seq, 1, length(seq))
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
