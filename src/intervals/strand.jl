# Strand
# ======
#
# Strand type for genomic ranges.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

bitstype 8 Strand

Base.convert(::Type{Strand}, strand::UInt8) = reinterpret(Strand, strand)
Base.convert(::Type{UInt8}, strand::Strand) = reinterpret(UInt8, strand)

const STRAND_NA   = convert(Strand, 0b000)
const STRAND_POS  = convert(Strand, 0b001)
const STRAND_NEG  = convert(Strand, 0b010)
const STRAND_BOTH = convert(Strand, 0b011)

function Base.show(io::IO, strand::Strand)
    if strand == STRAND_NA
        print(io, "?")
    elseif strand == STRAND_POS
        print(io, "+")
    elseif strand == STRAND_NEG
        print(io, "-")
    elseif strand == STRAND_BOTH
        print(io, ".")
    else
        print(io, "(undefined strand)")
    end
end

Base.isless(a::Strand, b::Strand) = convert(UInt8, a) < convert(UInt8, b)

function Base.convert(::Type{Strand}, strand::Char)
    if strand == '+'
        return STRAND_POS
    elseif strand == '-'
        return STRAND_NEG
    elseif strand == '.'
        return STRAND_BOTH
    elseif strand == '?'
        return STRAND_NA
    else
        error("$(strand) is not a valid strand")
    end
end
