# Constructors
# ============
#
# Constructor methods for Biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

function (::Type{BioSequence{A}}){A<:Alphabet}(len::Integer)
    return BioSequence{A}(Vector{UInt64}(seq_data_len(A, len)), 1:len, false)
end

BioSequence(::Type{DNA}) = DNASequence()
BioSequence(::Type{RNA}) = RNASequence()
BioSequence(::Type{AminoAcid}) = AminoAcidSequence()
BioSequence(::Type{Char}) = CharSequence()

function BioSequence()
    return BioSequence{VoidAlphabet}(Vector{UInt64}(0), 0:-1, false)
end

function (::Type{BioSequence{A}}){A<:Alphabet}(
        src::Union{AbstractString,AbstractVector},
        startpos::Integer=1,
        stoppos::Integer=length(src))
    len = stoppos - startpos + 1
    seq = BioSequence{A}(len)
    return encode_copy!(seq, 1, src, startpos, len)
end

# create a subsequence
function BioSequence{A,T<:Integer}(other::BioSequence{A}, part::UnitRange{T})
    checkbounds(other, part)
    start = other.part.start + part.start - 1
    stop = start + length(part) - 1
    subseq = BioSequence{A}(other.data, start:stop, true)
    other.shared = true
    return subseq
end

function (::Type{BioSequence{A}}){A}(other::BioSequence{A}, part::UnitRange)
    return BioSequence(other, part)
end

# concatenate chunks
function (::Type{BioSequence{A}}){A}(chunks::BioSequence{A}...)
    len = 0
    for chunk in chunks
        len += length(chunk)
    end
    seq = BioSequence{A}(len)
    offset = 1
    for chunk in chunks
        copy!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

function Base.repeat{A}(chunk::BioSequence{A}, n::Integer)
    seq = BioSequence{A}(length(chunk) * n)
    offset = 1
    for _ in 1:n
        copy!(seq, offset, chunk, 1)
        offset += length(chunk)
    end
    return seq
end

# operators for concat and repeat
Base.:*{A}(chunk::BioSequence{A}, chunks::BioSequence{A}...) =
    BioSequence{A}(chunk, chunks...)
Base.:^(chunk::BioSequence, n::Integer) = repeat(chunk, n)

function Base.similar{A}(seq::BioSequence{A}, len::Integer=length(seq))
    return BioSequence{A}(len)
end
