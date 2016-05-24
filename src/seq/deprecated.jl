import Base.@deprecate, Base.@deprecate_binding
import Base.depwarn

# DEPRECATED: defined for compatibility
type NucleotideSequence{T<:Nucleotide} end

@deprecate NucleotideSequence BioSequence
@deprecate NucleotideSequence(::Type{DNANucleotide}) DNASequence
@deprecate NucleotideSequence(::Type{RNANucleotide}) RNASequence
@compat begin
    (::Type{NucleotideSequence{DNANucleotide}})() = DNASequence()
    (::Type{NucleotideSequence{RNANucleotide}})() = RNASequence()
    (::Type{NucleotideSequence{DNANucleotide}})(
        seq::Union{AbstractVector{DNANucleotide},AbstractString}) = DNASequence(seq)
    (::Type{NucleotideSequence{RNANucleotide}})(
        seq::Union{AbstractVector{RNANucleotide},AbstractString}) = RNASequence(seq)
end
Base.convert(::Type{NucleotideSequence}, seq::DNAKmer) = DNASequence(seq)
Base.convert(::Type{NucleotideSequence}, seq::RNAKmer) = RNASequence(seq)
Base.convert(::Type{NucleotideSequence{DNANucleotide}}, seq::Union{AbstractVector{DNANucleotide},AbstractString,DNAKmer}) = DNASequence(seq)
Base.convert(::Type{NucleotideSequence{RNANucleotide}}, seq::Union{AbstractVector{RNANucleotide},AbstractString,RNAKmer}) = RNASequence(seq)

@deprecate NucleotideCounts(seq) Composition(seq)
@deprecate npositions(seq::BioSequence) ambiguous_positions(seq)
