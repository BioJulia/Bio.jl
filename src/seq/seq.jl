module Seq

using Compat
using Docile
using Base.Intrinsics

import Base: convert, getindex, show, length, start, next, done, copy, reverse,
             show, endof, ==, read, read!, isless
             show, endof, ==, isless

export Nucleotide, DNANucleotide, RNANucleotide,
       DNA_A, DNA_C, DNA_G, DNA_T, DNA_N,
       RNA_A, RNA_C, RNA_G, RNA_U, RNA_N,
       NucleotideSequence, DNASequence, RNASequence, @dna_str, @rna_str,
       complement, reverse_complement, mismatches, npositions, hasn, eachsubseq,
       eachkmer, NucleotideCounts, Kmer, DNAKmer, RNAKmer, dnakmer, rnakmer,
       kmer, AminoAcid, AminoAcidSequence, @aa_str, translate, ncbi_trans_table,
       AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L,
       AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_X,
       FASTA, FASTQ


abstract FileFormat
immutable FASTA <: FileFormat end
immutable FASTQ <: FileFormat end


abstract Sequence


include("nucleotide.jl")
include("aminoacid.jl")
include("geneticcode.jl")
include("alphabet.jl")
include("quality.jl")


# A sequence record is a named sequence with attached metadata.
@doc """
`SeqRecord{S, T}` is a type holding a named sequence of type `S`, along with
some arbitrary metadata of type `T`.
""" ->
type SeqRecord{S, T}
    name::String
    seq::S
    metadata::T

    function SeqRecord(name, seq::S, metadata::T)
        return new(name, seq, metadata)
    end

    function SeqRecord()
        return new("", S(), T())
    end
end


# Degelgate sequence operations
@doc """
Return a `SeqRecord` holding just the nucleotide at position `i`.
""" ->
function Base.getindex(seqrec::SeqRecord, i::Integer)
    return SeqRecord(seqreq.name, seqrec.seq[i], seqreq.metadata)
end

@doc """
Return a `SeqRecord` holding the specified subsequence
""" ->
function Base.getindex(seqrec::SeqRecord, r::UnitRange)
    return SeqRecord(seqrec.name, seqrec.seq[r], seqrec.metadata)
end


typealias DNASeqRecord{T}       SeqRecord{DNASequence, T}
typealias RNASeqRecord{T}       SeqRecord{RNASequence, T}
typealias AminoAcidSeqRecord{T} SeqRecord{RNASequence, T}


# Parsing of various file types
include("fasta.jl")
include("fastq.jl")

end # module Seq
