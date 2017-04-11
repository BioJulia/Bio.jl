# v0.4
# ----

include("hts_old/hts.jl")


# v0.3
# ----
@deprecate seq2ref(i::Integer, aln::Alignment) seq2ref(aln::Alignment, i::Integer)
@deprecate ref2seq(i::Integer, aln::Alignment) ref2seq(aln::Alignment, i::Integer)
@deprecate seq2ref(i::Integer, alnseq::AlignedSequence) seq2ref(aln::AlignedSequence, i::Integer)
@deprecate ref2seq(i::Integer, alnseq::AlignedSequence) ref2seq(aln::AlignedSequence, i::Integer)
@deprecate seq2ref(i::Integer, aln::PairwiseAlignment) seq2ref(aln::PairwiseAlignment, i::Integer)
@deprecate ref2seq(i::Integer, aln::PairwiseAlignment) ref2seq(aln::PairwiseAlignment, i::Integer)
