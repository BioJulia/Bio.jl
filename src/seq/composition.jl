# Composition
# ===========

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
    for x in seq
        comp.counts[encode(A, x)+1] += 1
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

function Base.getindex{A}(comp::Composition{A}, x)
    i = encode(A, x) + 1
    if i > endof(comp.counts)
        return 0
    else
        return Int(comp.counts[i])
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
