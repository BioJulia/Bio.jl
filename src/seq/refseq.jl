# ReferenceSequence
# =================
#
# DNA sequence for reference genomes.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
Reference Sequence

Reference sequence is a sequence of A/C/G/T/N. In the internals, it compresses
`N` positions and consumes less than three bits per base. Unlike `BioSequence`,
reference sequences are immutable and hence no modifyting operators are
provided.
"""
immutable ReferenceSequence <: Sequence
    data::Vector{UInt64}  # 2-bit encoding of A/C/G/T nucloetides
    nmask::NMask          # positions of N
    part::UnitRange{Int}  # interval within `data` defining the (sub)sequence
end

Base.length(seq::ReferenceSequence) = length(seq.part)
Base.eltype(::Type{ReferenceSequence}) = DNA
Base.summary(seq::ReferenceSequence) = string(length(seq), "nt Reference Sequence")

function Base.copy(seq::ReferenceSequence)
    return ReferenceSequence(copy(seq.data), copy(seq.nmask), seq.part)
end

function ReferenceSequence()
    return ReferenceSequence(UInt64[], NMask(), 1:0)
end

function ReferenceSequence{T<:Integer}(seq::ReferenceSequence, part::UnitRange{T})
    ReferenceSequence(seq.data, seq.nmask, part)
end

function ReferenceSequence(src::Vector{UInt8}, startpos::Integer=1,
                           len::Integer=length(src))
    return encode(src, startpos, len)
end

function Base.convert{A<:DNAAlphabet}(::Type{ReferenceSequence}, seq::BioSequence{A})
    data = Vector{UInt64}(cld(length(seq), 32))
    nmask = falses(length(seq))
    i = 1
    for j in 1:endof(data)
        x = UInt64(0)
        r = 0
        while r < 64 && i ≤ endof(seq)
            nt = seq[i]
            if nt == DNA_A
                x |= convert(UInt64, 0) << r
            elseif nt == DNA_C
                x |= convert(UInt64, 1) << r
            elseif nt == DNA_G
                x |= convert(UInt64, 2) << r
            elseif nt == DNA_T
                x |= convert(UInt64, 3) << r
            elseif nt == DNA_N
                nmask[i] = true
            else
                throw(ArgumentError("invalid symbol $(seq[i]) ∉ {A,C,G,T,N} at $i"))
            end
            i += 1
            r += 2
        end
        data[j] = x
    end
    return ReferenceSequence(data, NMask(nmask), 1:length(seq))
end

function Base.convert(::Type{ReferenceSequence}, str::AbstractString)
    if !isascii(str)
        throw(ArgumentError("attempt to convert a non-ASCII string to ReferenceSequence"))
    end
    return encode([UInt8(char) for char in str], 1, length(str))
end

function Base.convert(::Type{DNASequence}, seq::ReferenceSequence)
    bioseq = DNASequence(length(seq))
    for i in 1:endof(seq)
        bioseq[i] = seq[i]
    end
    return bioseq
end

function Base.convert{S<:AbstractString}(::Type{S}, seq::ReferenceSequence)
    return S([Char(nt) for nt in seq])
end

function bitindex(seq::ReferenceSequence, i::Integer)
    return BitIndex((i + first(seq.part) - 2) << 1)
end

# create ReferenceSequence object from the ascii-encoded `data`
function encode(src::Vector{UInt8}, from::Integer, len::Integer)
    data = zeros(UInt64, cld(len, 32))
    nmask = falses(len)
    next = BitIndex(1, 2)
    stop = BitIndex(len + 1, 2)
    i = from
    while next < stop
        x = UInt64(0)
        j = index(next)
        while index(next) == j && next < stop
            # FIXME: Hotspot
            char = convert(Char, src[i])
            nt = convert(DNA, char)
            if !isambiguous(nt)
                x |= UInt64(encode(DNAAlphabet{2}, nt)) << offset(next)
            elseif nt == DNA_N
                nmask[i] = true
            else
                error("'", char, "'", " is not allowed")
            end
            i += 1
            next += 2
        end
        data[j] = x
    end
    return ReferenceSequence(data, nmask, 1:len)
end

function Base.checkbounds(seq::ReferenceSequence, part::UnitRange)
    if isempty(part) || (1 ≤ first(part) && last(part) ≤ endof(seq))
        return true
    end
    throw(BoundsError(seq, part))
end

@inline function inbounds_getindex(seq::ReferenceSequence, i::Integer)
    if seq.nmask[i + first(seq.part) - 1]
        return DNA_N
    else
        j = bitindex(seq, i)
        return DNA(0x01 << ((seq.data[index(j)] >> offset(j)) & 0b11))
    end
end

function Base.getindex{T<:Integer}(seq::ReferenceSequence, part::UnitRange{T})
    checkbounds(seq, part)
    return ReferenceSequence(seq, part)
end

function find_next_ambiguous(seq::ReferenceSequence, i::Integer)
    return findnextn(seq.nmask, i)
end
