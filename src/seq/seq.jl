
module Seq

using Base.Intrinsics
using IntervalTrees
import Base: convert, getindex, show, length, start, next, done

export NucleotideSequence, DNANucleotide, RNANucleotide,
       DNASequence, RNASequence, @dna_str, @rna_str


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


const DNA_A = convert(DNANucleotide, 0b000)
const DNA_T = convert(DNANucleotide, 0b001)
const DNA_C = convert(DNANucleotide, 0b010)
const DNA_G = convert(DNANucleotide, 0b011)
const DNA_N = convert(DNANucleotide, 0b100)

const RNA_A = convert(RNANucleotide, 0b000)
const RNA_U = convert(RNANucleotide, 0b001)
const RNA_C = convert(RNANucleotide, 0b010)
const RNA_G = convert(RNANucleotide, 0b011)
const RNA_N = convert(RNANucleotide, 0b100)


const dna_to_char = ['A', 'T', 'C', 'G', 'N']

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


const rna_to_char = ['A', 'U', 'C', 'G', 'N']

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

type NucleotideSequence{T <: Nucleotide}
    # 2-bit encoded sequence
    data::Vector{Uint64}

    # length of the sequence
    len::Int

    # 'N' intervals
    ns::IntervalTree{Int, Bool}

    function NucleotideSequence{T}()
        return new(zeros(Uint64, 0), len)
    end

    function NucleotideSequence(seq::String)
        data = zeros(Uint64, div(length(seq), 32) + 1)
        ns = IntervalTree{Int, Bool}()

        nstart = 0
        for (i, nt) in enumerate(seq)
            if nstart != 0
                if nt != 'N' && nt != 'n'
                    ns[(nstart, i-1)] = true
                    nstart = 0
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

        return new(data, length(seq), ns)
    end

    # TODO: A constructor that looks for Ts or Us to determine wether
    # the sequence is DNA or RNA.
end


typealias DNASequence NucleotideSequence{DNANucleotide}
typealias RNASequence NucleotideSequence{RNANucleotide}


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
    return seq.len
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
    if i > seq.len || i < 1
        error(BoundsError())
    end
    d, r = divrem(i - 1, 32)
    return convert(T, convert(Uint8, (seq.data[d + 1] >> (2*r)) & 0b11))
end


function convert(::Type{DNASequence}, seq::String)
    return DNASequence(seq)
end


function convert(::Type{RNASequence}, seq::String)
    return DNASequence(seq)
end


function convert(::Type{String}, seq::NucleotideSequence)
    # TODO
end


end
