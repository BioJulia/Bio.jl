# Composition
# ===========
#
# Sequence composition counter.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

type Composition{A<:Alphabet}
    counts::Vector{Int}
    function Composition()
        new(zeros(Int, length(alphabet(A))))
    end
end

typealias DNANucleotideCounts{n} Composition{DNAAlphabet{n}}
typealias RNANucleotideCounts{n} Composition{RNAAlphabet{n}}
typealias AminoAcidCounts Composition{AminoAcidAlphabet}

function Composition{A}(seq::BioSequence{A})
    comp = Composition{A}()
    @inbounds for x in seq
        comp.counts[encode(A, x)+1] += 1
    end
    return comp
end

function Composition(seq::ReferenceSequence)
    comp = Composition{DNAAlphabet{4}}()
    @inbounds for x in seq
        comp.counts[encode(DNAAlphabet{4}, x)+1] += 1
    end
    return comp
end

function Composition{K}(kmer::DNAKmer{K})
    A = DNAAlphabet{2}
    comp = Composition{A}()
    x = UInt64(kmer)
    @inbounds begin
        comp.counts[encode(A, DNA_A)+1] = count_a(x) - (32 - K)
        comp.counts[encode(A, DNA_C)+1] = count_c(x)
        comp.counts[encode(A, DNA_G)+1] = count_g(x)
        comp.counts[encode(A, DNA_T)+1] = count_t(x)
    end
    return comp
end

function Composition{K}(kmer::RNAKmer{K})
    A = RNAAlphabet{2}
    comp = Composition{A}()
    x = UInt64(kmer)
    @inbounds begin
        comp.counts[encode(A, RNA_A)+1] = count_a(x) - (32 - K)
        comp.counts[encode(A, RNA_C)+1] = count_c(x)
        comp.counts[encode(A, RNA_G)+1] = count_g(x)
        comp.counts[encode(A, RNA_U)+1] = count_t(x)
    end
    return comp
end

"""
    composition(seq)

Calculate composition of biological symbols in `seq`.
"""
function composition(seq::Sequence)
    return Composition(seq)
end

function Base.getindex{A}(comp::Composition{A}, x)
    try
        i = encode(A, x) + 1
        return Int(comp.counts[i])
    catch ex
        if isa(ex, EncodeError{A})
            return 0
        end
        rethrow()
    end
end

Base.summary(::DNANucleotideCounts) = "DNA Nucleotide Counts"
Base.summary(::RNANucleotideCounts) = "RNA Nucleotide Counts"
Base.summary(::AminoAcidCounts)     = "Amino Acid Counts"

function Base.show{A}(io::IO, comp::Composition{A})
    print(io, summary(comp), ':')
    for x in alphabet(A)
        print(io, "\n  ", x, " => ", comp[x])
    end
end
