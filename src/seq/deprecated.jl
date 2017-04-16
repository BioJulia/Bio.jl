import Base.@deprecate, Base.@deprecate_binding
import Base.depwarn

# v0.4
# ----

@deprecate_binding Nucleotide NucleicAcid
@deprecate_binding DNANucleotide DNA
@deprecate_binding RNANucleotide RNA

include("fasta_old/fasta.jl")
include("fastq_old/fastq.jl")
include("twobit_old/twobit.jl")

export
    FASTAReader,
    FASTAWriter,
    FASTASeqRecord,
    FASTQReader,
    FASTQWriter,
    FASTQSeqRecord,
    TwoBitReader,
    TwoBitWriter


# v0.3
# ----

@deprecate ispalindrome ispalindromic

"Range type of biological symbols (especially `DNA`, `RNA`, and `AminoAcid`)."
immutable SymbolRange{T} <: Range{T}
    start::UInt8
    stop::UInt8

    function SymbolRange(start::T, stop::T)
        warn("SymbolRange is deprecated")
        if start > stop
            # normalize empty range
            start = convert(T, 0x01)
            stop = convert(T, 0x00)
        end
        return new(start, stop)
    end
end

function SymbolRange{T}(start::T, stop::T)
    return SymbolRange{T}(start, stop)
end

function Base.show(io::IO, r::SymbolRange)
    show(io, first(r))
    write(io, ':')
    show(io, last(r))
end


# Iterator
# --------

Base.start(r::SymbolRange) = r.start
Base.done(r::SymbolRange, i) = i > r.stop
Base.next(r::SymbolRange, i) = reinterpret(eltype(r), i), i + 0x01

Base.size(r::SymbolRange) = (r.stop - r.start + 1,)
Base.first(r::SymbolRange) = convert(eltype(r), r.start)
Base.last(r::SymbolRange) = convert(eltype(r), r.stop)
Base.in{T}(x::T, r::SymbolRange{T}) = first(r) ≤ x ≤ last(r)


# Indexing
# --------

function Base.getindex(r::SymbolRange, i::Integer)
    checkbounds(r, i)
    return convert(eltype(r), r.start + i - 1)
end

function Base.getindex{T}(r::SymbolRange{T}, ir::UnitRange)
    checkbounds(r, ir)
    if isempty(ir)
        return SymbolRange(T(0x01), T(0x00))
    else
        return SymbolRange(r[first(ir)], r[last(ir)])
    end
end
