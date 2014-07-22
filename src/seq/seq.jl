
module Seq

using Base.Intrinsics

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof

export Nucleotide, DNANucleotide, RNANucleotide,
       DNA_A, DNA_C, DNA_G, DNA_T, DNA_N,
       RNA_A, RNA_C, RNA_G, RNA_T, RNA_N,
       NucleotideSequence, DNASequence, RNASequence, @dna_str, @rna_str,
       complement, reverse_complement, mismatches, npositions, eachsubseq,
       eachkmer, nucleotide_count, Kmer, DNAKmer, RNAKmer, dnakmer, rnakmer


include("nucleotide.jl")
include("kmer.jl")
# TODO: aminoacid
# TODO: pattern (pseudo-sequences with ambiguity codes, etc)

end # module Seq

