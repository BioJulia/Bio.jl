# Bio.Seq
# =======
#
# Module for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

module Seq

export
    NucleicAcid,
    DNA,
    RNA,
    DNA_A,
    DNA_C,
    DNA_G,
    DNA_T,
    DNA_M,
    DNA_R,
    DNA_W,
    DNA_S,
    DNA_Y,
    DNA_K,
    DNA_V,
    DNA_H,
    DNA_D,
    DNA_B,
    DNA_N,
    DNA_Gap,
    ACGT,
    RNA_A,
    RNA_C,
    RNA_G,
    RNA_U,
    RNA_M,
    RNA_R,
    RNA_W,
    RNA_S,
    RNA_Y,
    RNA_K,
    RNA_V,
    RNA_H,
    RNA_D,
    RNA_B,
    RNA_N,
    RNA_Gap,
    ACGU,
    iscompatible,
    isambiguous,
    iscertain,
    isgap,
    ispurine,
    ispyrimidine,
    Sequence,
    BioSequence,
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    CharSequence,
    NucleicAcidSequence,
    SeqRecord,
    seqname,
    sequence,
    metadata,
    @dna_str,
    @rna_str,
    @aa_str,
    @char_str,
    @biore_str,
    @prosite_str,
    matched,
    captured,
    alphabet,
    gap,
    complement,
    complement!,
    reverse_complement,
    reverse_complement!,
    mismatches,
    ispalindromic,
    hasambiguity,
    isrepetitive,
    ambiguous_positions,
    gc_content,
    SequenceGenerator,
    randdnaseq,
    randrnaseq,
    randaaseq,
    canonical,
    neighbors,
    eachkmer,
    each,
    Composition,
    composition,
    NucleicAcidCounts,
    Kmer,
    DNAKmer,
    RNAKmer,
    DNACodon,
    RNACodon,
    translate,
    ncbi_trans_table,
    AminoAcid,
    AA_A,
    AA_R,
    AA_N,
    AA_D,
    AA_C,
    AA_Q,
    AA_E,
    AA_G,
    AA_H,
    AA_I,
    AA_L,
    AA_K,
    AA_M,
    AA_F,
    AA_P,
    AA_S,
    AA_T,
    AA_W,
    AA_Y,
    AA_V,
    AA_O,
    AA_U,
    AA_B,
    AA_J,
    AA_Z,
    AA_X,
    AA_Term,
    AA_Gap,
    FASTA,
    FASTAReader,
    FASTAWriter,
    FASTASeqRecord,
    FASTQReader,
    FASTQWriter,
    FASTQSeqRecord,
    TwoBitReader,
    TwoBitWriter,
    Alphabet,
    DNAAlphabet,
    RNAAlphabet,
    NucleicAcidAlphabets,
    AminoAcidAlphabet,
    CharAlphabet,
    NucleicAcidAlphabet,
    ExactSearchQuery,
    ApproximateSearchQuery,
    approxsearch,
    approxsearchindex,
    approxrsearch,
    approxrsearchindex,
    ReferenceSequence,
    Demultiplexer,
    demultiplex,
    seqmatrix,
    majorityvote,
    tryread!

import Automa
import Automa.RegExp: @re_str
import Bio.Ragel: Ragel, tryread!
import Bio.StringFields: StringField
import BufferedStreams: BufferedStreams, BufferedInputStream, BufferedOutputStream
import Combinatorics
import IndexableBitVectors
import Iterators
importall Bio

"""
    alphabet(typ)

Return an iterator of symbols of `typ`.

`typ` is one of `DNA`, `RNA`, or `AminoAcid`.
"""
function alphabet end

"""
    gap(typ)

Return the gap symbol of `typ`.

`typ` is one of `DNA`, `RNA`, `AminoAcid`, or `Char`.
"""
function gap end

gap(::Type{Char}) = '-'

include("nucleicacid.jl")
include("aminoacid.jl")
include("alphabet.jl")
include("bitindex.jl")
include("sequence.jl")
include("bioseq/bioseq.jl")
include("bioseq/hash.jl")
include("bioseq/randseq.jl")
include("kmer.jl")
include("nmask.jl")
include("refseq.jl")
include("eachkmer.jl")
include("composition.jl")
include("geneticcode.jl")
include("seqrecord.jl")
include("demultiplexer.jl")

# Parsing file types
include("fasta/fasta.jl")
include("fastq/fastq.jl")
include("twobit/twobit.jl")

include("search/exact.jl")
include("search/approx.jl")
include("search/re.jl")

include("deprecated.jl")

end  # module Bio.Seq
