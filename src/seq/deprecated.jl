import Base.@deprecate, Base.@deprecate_binding
import Base.depwarn

# v0.2
# ----

"""
    DNAKmerCounts{k}(seq::Sequence[, step=1])
    RNAKmerCounts{k}(seq::Sequence[, step=1])

Count ocurrences of short (≤ 32) `k`-mers in a sequence.

This method uses a dense table to store counts, requiring O(4^k) memory, so is
not recommended for large k-mer sizes.

# Arguments
* `seq`: A nucleotide sequence.
* `step=1`: K-mers counted are separated by this many nucleotides.
"""
immutable KmerCounts{T,K}
    data::Vector{UInt32}

    function KmerCounts(seq::Sequence, step::Integer=1)
        data = zeros(UInt32, 4^K)
        @inbounds for (_, x) in each(Kmer{T,K}, seq, step)
            data[convert(UInt64, x) + 1] += 1
        end
        return new(data)
    end
end

typealias DNAKmerCounts{K} KmerCounts{DNANucleotide,K}
typealias RNAKmerCounts{K} KmerCounts{RNANucleotide,K}

function Base.getindex{T,K}(counts::KmerCounts{T,K}, x::Kmer{T,K})
    @inbounds c = counts.data[convert(UInt64, x) + 1]
    return c
end

function Base.show{T,K}(io::IO, counts::KmerCounts{T,K})
    println(io, (T == DNANucleotide ? "DNA" : "RNA"), "KmerCounts{", K, "}:")
    for x in UInt64(1):UInt64(4^K)
        s = string(Kmer{T,K}(x - 1))
        println(io, "  ", s, " => ", counts.data[x])
    end
end

@deprecate KmerCounts composition


# v0.3
# ----

@deprecate ispalindrome ispalindromic

"""
Range type of biological symbols (especially `DNANucleotide`, `RNANucleotide`,
and `AminoAcid`).
"""
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


immutable FASTA <: Bio.IO.FileFormat end
immutable FASTQ <: Bio.IO.FileFormat end
immutable TwoBit <: Bio.IO.FileFormat end

function Base.open(filepath::AbstractString, ::Type{FASTA})
    input = BufferedInputStream(open(filepath))
    indexpath = filepath * ".fai"
    if isfile(indexpath)
        return FASTAReader{BioSequence}(input, FASTAIndex(indexpath))
    else
        return FASTAReader{BioSequence}(input)
    end
end

function Base.open(filepath::AbstractString, mode::AbstractString, ::Type{FASTA};
                   width::Integer=60)
    io = open(filepath, mode)
    if mode[1] == 'r'
        return open(BufferedInputStream(io), FASTA)
    elseif mode[1] ∈ ('w', 'a')
        return FASTAWriter(io, width)
    end
    error("invalid open mode")
end

function Base.open{S}(input::BufferedInputStream, ::Type{FASTA},
                      ::Type{S}=BioSequence)
    return FASTAReader{S}(input)
end
