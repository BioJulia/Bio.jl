module Seq

export
    Nucleotide,
    DNANucleotide,
    RNANucleotide,
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
    Sequence,
    BioSequence,
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    CharSequence,
    NucleotideSequence,
    @dna_str,
    @rna_str,
    @aa_str,
    @char_str,
    alphabet,
    gap,
    complement,
    complement!,
    reverse_complement,
    reverse_complement!,
    mismatches,
    ambiguous_positions,
    npositions,
    hasn,
    canonical,
    neighbors,
    eachkmer,
    each,
    Composition,
    composition,
    NucleotideCounts,
    Kmer,
    DNAKmer,
    RNAKmer,
    dnakmer,
    rnakmer,
    kmer,
    KmerCounts,
    DNAKmerCounts,
    RNAKmerCounts,
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
    FASTQ,
    Alphabet,
    DNAAlphabet,
    RNAAlphabet,
    AminoAcidAlphabet,
    CharAlphabet,
    approxsearch,
    approxsearchindex,
    approxrsearch,
    approxrsearchindex

using BufferedStreams,
    Base.Intrinsics,
    Bio.StringFields,
    Bio.Ragel

using Bio:
    FileFormat,
    AbstractParser


abstract Sequence

# This is useful for obscure reasons. We use SeqRecord{Sequence} for reading
# sequence in an undetermined alphabet, but a consequence that we need to be
# able to construct a `Sequence`.
function Sequence()
    return DNASequence()
end

"""
    alphabet(typ)

Return an iterator of symbols of `typ`.

`typ` is one of `DNANucleotide`, `RNANucleotide`, or `AminoAcid`.
"""
function alphabet end

"""
    gap(typ)

Return the gap symbol of `typ`.

`typ` is one of `DNANucleotide`, `RNANucleotide`, `AminoAcid`, or `Char`.
"""
function gap end

gap(::Type{Char}) = '-'

include("nucleotide.jl")
include("aminoacid.jl")
include("alphabet.jl")
include("bioseq.jl")
include("hash.jl")
include("kmer.jl")
include("eachkmer.jl")
include("kmercounts.jl")
include("composition.jl")
include("geneticcode.jl")
include("util.jl")
include("quality.jl")
include("seqrecord.jl")

# Parsing file types
include("naivebayes.jl")
include("fasta.jl")
include("fastq.jl")

include("search/exact.jl")
include("search/approx.jl")
include("search/search.jl")

# TODO: any way to write this elegantly?
import .RE:
    Regex,
    RegexMatch,
    matched,
    captured,
    @dnapat_str,
    @rnapat_str,
    @aapat_str,
    @prosite_str
export
    matched,
    captured,
    @dnapat_str,
    @rnapat_str,
    @aapat_str,
    @prosite_str

>>>>>>> add moved regular expression

# DEPRECATED: defined just for compatibility
type NucleotideSequence{T<:Nucleotide} end

Base.call(::Type{NucleotideSequence{DNANucleotide}}) = DNASequence()
Base.call(::Type{NucleotideSequence{RNANucleotide}}) = RNASequence()
Base.call(::Type{NucleotideSequence}, ::Type{DNANucleotide}) = DNASequence()
Base.call(::Type{NucleotideSequence}, ::Type{RNANucleotide}) = RNASequence()
Base.call(::Type{NucleotideSequence{DNANucleotide}}, seq::Union{AbstractVector{DNANucleotide},AbstractString}) = DNASequence(seq)
Base.call(::Type{NucleotideSequence{RNANucleotide}}, seq::Union{AbstractVector{RNANucleotide},AbstractString}) = RNASequence(seq)

Base.convert(::Type{NucleotideSequence}, seq::DNAKmer) = DNASequence(seq)
Base.convert(::Type{NucleotideSequence}, seq::RNAKmer) = RNASequence(seq)
Base.convert(::Type{NucleotideSequence{DNANucleotide}}, seq::Union{AbstractVector{DNANucleotide},AbstractString,DNAKmer}) = DNASequence(seq)
Base.convert(::Type{NucleotideSequence{RNANucleotide}}, seq::Union{AbstractVector{RNANucleotide},AbstractString,RNAKmer}) = RNASequence(seq)

# DEPRECATED
NucleotideCounts(seq) = Composition(seq)

# DEPRECATED: use ambiguous_positions
npositions(seq::BioSequence) = ambiguous_positions(seq)

end # module Seq
