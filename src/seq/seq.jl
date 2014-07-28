
module Seq

using Base.Intrinsics

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof

export Nucleotide, DNANucleotide, RNANucleotide,
       DNA_A, DNA_C, DNA_G, DNA_T, DNA_N,
       RNA_A, RNA_C, RNA_G, RNA_T, RNA_N,
       NucleotideSequence, DNASequence, RNASequence, @dna_str, @rna_str,
       complement, reverse_complement, mismatches, npositions, hasn, eachsubseq,
       eachkmer, nucleotide_count, Kmer, DNAKmer, RNAKmer, dnakmer, rnakmer,
       AminoAcid, AminoAcidSequence, @aa_str, translate, ncbi_trans_table,
       AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L,
       AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_X


include("nucleotide.jl")
include("kmer.jl")
include("aminoacid.jl")
# TODO: pattern (pseudo-sequences with ambiguity codes, etc)

end # module Seq

