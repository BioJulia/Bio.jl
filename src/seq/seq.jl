
module Seq

using Base.Intrinsics
using IntervalTrees

import Base: convert, getindex, show, length, start, next, done, copy

export Nucleotide, DNANucleotide, RNANucleotide,
       NucleotideSequence, DNASequence, RNASequence, @dna_str, @rna_str


# Nucleotides
# -----------


abstract Nucleotide
bitstype 8 DNANucleotide <: Nucleotide
bitstype 8 RNANucleotide <: Nucleotide


function convert(::Type{DNANucleotide}, nt::Uint8)
    return box(DNANucleotide, unbox(Uint8, nt))
end


function convert(::Type{RNANucleotide}, nt::Uint8)
    return box(RNANucleotide, unbox(Uint8, nt))
end


function convert(::Type{Uint8}, nt::DNANucleotide)
    return box(Uint8, unbox(DNANucleotide, nt))
end


function convert(::Type{Uint8}, nt::RNANucleotide)
    return box(Uint8, unbox(RNANucleotide, nt))
end


function convert{T <: Integer}(::Type{T}, nt::DNANucleotide)
    return convert(T, convert(Uint8, nt))
end


function convert{T <: Integer}(::Type{T}, nt::RNANucleotide)
    return convert(T, convert(Uint8, nt))
end


# TODO: Should we use ape-style bit masks for these?
const DNA_A = convert(DNANucleotide, 0b000)
const DNA_C = convert(DNANucleotide, 0b001)
const DNA_G = convert(DNANucleotide, 0b010)
const DNA_T = convert(DNANucleotide, 0b011)
const DNA_N = convert(DNANucleotide, 0b100)

const RNA_A = convert(RNANucleotide, 0b000)
const RNA_C = convert(RNANucleotide, 0b001)
const RNA_G = convert(RNANucleotide, 0b010)
const RNA_U = convert(RNANucleotide, 0b011)
const RNA_N = convert(RNANucleotide, 0b100)


const dna_to_char = ['A', 'C', 'G', 'T', 'N']

function convert(::Type{Char}, nt::DNANucleotide)
    return dna_to_char[convert(Uint8, nt) + 1]
end

# TODO: it may be faster to look things up in an array
function convert(::Type{DNANucleotide}, nt::Char)
    if nt == 'A' || nt == 'a'
        return DNA_A
    elseif nt == 'C' || nt == 'c'
        return DNA_C
    elseif nt == 'T' || nt == 't'
        return DNA_T
    elseif nt == 'G' || nt == 'g'
        return DNA_G
    elseif nt == 'N' || nt == 'n'
        return DNA_N
    else
        error("$(nt) is not a valid DNA nucleotide")
    end
end


const rna_to_char = ['A', 'C', 'G', 'U', 'N']

function convert(::Type{Char}, nt::RNANucleotide)
    return rna_to_char[convert(Uint8, nt) + 1]
end


function convert(::Type{RNANucleotide}, nt::Char)
    if nt == 'A' || nt == 'a'
        return RNA_A
    elseif nt == 'C' || nt == 'c'
        return RNA_C
    elseif nt == 'U' || nt == 'U'
        return RNA_U
    elseif nt == 'G' || nt == 'g'
        return RNA_G
    elseif nt == 'N' || nt == 'n'
        return RNA_N
    else
        error("$(nt) is not a valid RNA nucleotide")
    end
end


function show(io::IO, nt::DNANucleotide)
    write(io, convert(Char, nt))
end


function show(io::IO, nt::RNANucleotide)
    write(io, convert(Char, nt))
end


# Nucleotide Sequences
# --------------------

# How many Uint64s are needed to represent a sequence of length n
function seq_data_len(n::Integer)
    d, r = divrem(n, 32)
    return d + (r > 0 ? 1 : 0)
end

type NucleotideSequence{T <: Nucleotide}
    # 2-bit encoded sequence
    data::Vector{Uint64}

    # 'N' intervals
    ns::IntervalTree{Int, Bool}

    # interval within data defining the (sub)sequence
    part::UnitRange{Int}

    # Construct from raw components
    function NucleotideSequence(data::Vector{Uint64}, ns::IntervalTree{Int, Bool},
                                part::UnitRange)
        return new(data, ns, part)
    end

    # Construct a subsequence of another nucleotide sequence
    function NucleotideSequence(other::NucleotideSequence, part::UnitRange)
        start = other.part.start + part.start - 1
        stop = start + length(part) - 1
        @assert start >= other.part.start
        @assert stop <= other.part.stop
        return new(other.data, other.ns, part)
    end

    # Construct an empty sequence
    function NucleotideSequence()
        return new(zeros(Uint64, 0), IntervalTree{Int, Bool}(), 1:0)
    end

    # Construct a sequence from a string
    function NucleotideSequence(seq::String)
        data = zeros(Uint64, seq_data_len(length(seq)))
        ns = IntervalTree{Int, Bool}()

        nstart = 0
        for (i, nt) in enumerate(seq)
            if nstart != 0
                if nt != 'N' && nt != 'n'
                    ns[(nstart, i-1)] = true
                    nstart = 0
                else
                    continue
                end
            elseif nt == 'N' || nt == 'n'
                nstart = i
                continue
            end

            d, r = divrem(i - 1, 32)
            data[d + 1] |= convert(Uint64, convert(T, nt)) << (2*r)
        end

        if nstart != 0
            ns[(nstart, length(seq))] = true
        end

        return new(data, ns, 1:length(seq))
    end

    # TODO: A constructor that looks for Ts or Us to determine wether
    # the sequence is DNA or RNA.
end


typealias DNASequence NucleotideSequence{DNANucleotide}
typealias RNASequence NucleotideSequence{RNANucleotide}


function copy{T}(seq::NucleotideSequence{T})
    data = zeros(Uint64, seq_data_len(length(seq)))
    d0, r0 = divrem(seq.part.start - 1, 32)

    h = 64 - 2*r0
    k = 2*r0

    j = d0
    for i in 1:length(data)
        data[i] |= seq.data[j] >>> k

        j += 1
        if j > length(seq.data)
            break
        end

        data[i] |= seq.data[j] << h
    end

    ns = IntervalTree{Int, Bool}()
    for (a, b) in intersect(seq.ns, (seq.part.start, seq.part.stop))
        ns[(a,b)] = true
    end

    return NucleotideSequence{T}(data, ns, 1:length(seq))
end

# Iterating throug nucleotide sequences
function start(seq::NucleotideSequence)
    return 1
end


function next(seq::DNASequence, i)
    value = hasintersection(seq.ns, i) ? DNA_N : seq[i];
    return value, i + 1
end


function next(seq::RNASequence, i)
    value = hasintersection(seq.ns, i) ? RNA_N : seq[i];
    return value, i + 1
end


function done(seq::NucleotideSequence, i)
    return i > length(seq)
end


# String decorator syntax to enable building sequence literals like:
#     dna"ACGTACGT" and rna"ACGUACGU"
macro dna_str(seq, flags...)
    return DNASequence(seq)
end

macro rna_str(seq, flags...)
    return RNASequence(seq)
end


function length(seq::NucleotideSequence)
    return length(seq.part)
end


function show{T}(io::IO, seq::NucleotideSequence{T})
    if T == DNANucleotide
        write(io, "dna\"")
    elseif T == RNANucleotide
        write(io, "rna\"")
    end

    len = length(seq)
    for nt in seq
        write(io, convert(Char, nt))
    end

    write(io, "\"")
end


function getindex{T}(seq::NucleotideSequence{T}, i::Integer)
    if i > length(seq) || i < 1
        error(BoundsError())
    end
    i += seq.part.start - 1
    if hasintersection(seq.ns, i)
        return convert(T, 0b0100) # bit representation for N
    else
        d, r = divrem(i - 1, 32)
        return convert(T, convert(Uint8, (seq.data[d + 1] >>> (2*r)) & 0b11))
    end
end


function getindex{T}(seq::NucleotideSequence{T}, r::UnitRange)
    return NucleotideSequence{T}(seq, r)
end


function convert(::Type{DNASequence}, seq::String)
    return DNASequence(seq)
end


function convert(::Type{RNASequence}, seq::String)
    return DNASequence(seq)
end


function convert(::Type{String}, seq::NucleotideSequence)
    return convert(String, [convert(Char, x) for x in seq])
end


# Transformations
# ---------------

# If seq is a sub-sequence do we want to take the complement of the full
# sequence?
function complement{T}(seq::NucleotideSequence{T})
    data = Array(Uint64, length(seq.data))


end



end
