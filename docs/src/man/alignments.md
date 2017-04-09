# Align: Sequence Alignments

```@meta
CurrentModule = Bio.Align
DocTestSetup = quote
    using Bio.Seq
    using Bio.Align
end
```

The `Align` module contains tools for computing and working with sequence
alignments.


## Representing alignments

The `Alignment` type can represent a wide variety of global or local sequence
alignments while facilitating efficient coordinate transformation.  Alignment
are always relative to a possibly unspecified reference sequence and represent a
series of [edit operations](https://en.wikipedia.org/wiki/Edit_distance)
performed on that reference to transform it to the query sequence.

To represent an alignment we use a series of "anchors" stored in the
`AlignmentAnchor` type. Anchors are form of run-length encoding alignment
operations, but rather than store an operation along with a length, we store the
end-point of that operation in both reference and query coordinates.

```julia
immutable AlignmentAnchor
    seqpos::Int
    refpos::Int
    op::Operation
end
```

Every alignment starts with a special `OP_START` operation which is used to give
the position in the reference and query prior to the start of the alignment, or
0, if the alignment starts at position 1.

For example, consider the following alignment:

                  0   4        9  12 15     19
                  |   |        |  |  |      |
        query:     TGGC----ATCATTTAACG---CAAG
    reference: AGGGTGGCATTTATCAG---ACGTTTCGAGAC
                  |   |   |    |     |  |   |
                  4   8   12   17    20 23  27

Using anchors we would represent this as the following series of anchors:
```julia
[
    AlignmentAnchor( 0,  4, OP_START),
    AlignmentAnchor( 4,  8, OP_MATCH),
    AlignmentAnchor( 4, 12, OP_DELETE),
    AlignmentAnchor( 9, 17, OP_MATCH),
    AlignmentAnchor(12, 17, OP_INSERT),
    AlignmentAnchor(15, 20, OP_MATCH),
    AlignmentAnchor(15, 23, OP_DELETE),
    AlignmentAnchor(19, 27, OP_MATCH)
]
```

An `Alignment` object can be created from a series of anchors:
```jlcon
julia> Alignment([
           AlignmentAnchor(0, 4, OP_START),
           AlignmentAnchor(4, 8, OP_MATCH),
           AlignmentAnchor(4, 12, OP_DELETE)
       ])
········
····----
```


### Operations

Alignment operations follow closely from those used in the [SAM/BAM
format](https://samtools.github.io/hts-specs/SAMv1.pdf) and are stored in the
`Operation` bitstype.

| Operation            | Operation Type     | Description                                                                     |
| :------------------- | :----------------- | :------------------------------------------------------------------------------ |
| `OP_MATCH`           | match              | non-specific match                                                              |
| `OP_INSERT`          | insert             | insertion into reference sequence                                               |
| `OP_DELETE`          | delete             | deletion from reference sequence                                                |
| `OP_SKIP`            | delete             | (typically long) deletion from the reference, e.g. due to RNA splicing          |
| `OP_SOFT_CLIP`       | insert             | sequence removed from the beginning or end of the query sequence but stored     |
| `OP_HARD_CLIP`       | insert             | sequence removed from the beginning or end of the query sequence and not stored |
| `OP_PAD`             | special            | not currently supported, but present for SAM/BAM compatibility                  |
| `OP_SEQ_MATCH`       | match              | match operation with matching sequence positions                                |
| `OP_SEQ_MISMATCH`    | match              | match operation with mismatching sequence positions                             |
| `OP_BACK`            | special            | not currently supported, but present for SAM/BAM compatibility                  |
| `OP_START`           | special            | indicate the start of an alignment within the reference and query sequence      |


## Aligned sequence

A sequence aligned to another sequence is represented by the `AlignedSequence`
type, which is a pair of the aligned sequence and an `Alignment` object.

The following example creates an aligned sequence object from a sequence and an
alignment:
```jlcon
julia> AlignedSequence(  # pass an Alignment object
           dna"ACGTAT",
           Alignment([
               AlignmentAnchor(0, 0, OP_START),
               AlignmentAnchor(3, 3, OP_MATCH),
               AlignmentAnchor(6, 3, OP_INSERT)
           ])
       )
···---
ACGTAT

julia> AlignedSequence(  # or pass a vector of anchors
           dna"ACGTAT",
           [
               AlignmentAnchor(0, 0, OP_START),
               AlignmentAnchor(3, 3, OP_MATCH),
               AlignmentAnchor(6, 3, OP_INSERT)
           ]
       )
···---
ACGTAT

```

If you already have an aligned sequence with gap symbols, it can be converted to
an `AlignedSequence` object by passing a reference sequence with it:
```jlcon
julia> seq = dna"ACGT--AAT--"
11nt DNA Sequence:
ACGT--AAT--

julia> ref = dna"ACGTTTAT-GG"
11nt DNA Sequence:
ACGTTTAT-GG

julia> AlignedSequence(seq, ref)
········-··
ACGT--AAT--

```


## Operating on alignments

```@docs
first
last
seq2ref
ref2seq
cigar
```


## Pairwise alignment

Pairwise alignment is a sequence alignment between two sequences.  The
`Bio.Align` module implements several pairwise alignment methods that maximize
alignment score or minimize alignment cost.

A pair of optimization type and score (or cost) model determines the best
alignments between two sequences. The next example uses a pair of
`GlobalAlignment` and `AffineGapScoreModel` to obtain the best alignment:

```jlcon
julia> problem = GlobalAlignment()
Bio.Align.GlobalAlignment()

julia> scoremodel = AffineGapScoreModel(
                  match=5,
                  mismatch=-4,
                  gap_open=-4,
                  gap_extend=-1
              )
Bio.Align.AffineGapScoreModel{Int64}:
       match = 5
    mismatch = -4
    gap_open = -4
  gap_extend = -1

julia> pairalign(problem, dna"CGGATTA", dna"GGTTTAC", scoremodel)
Bio.Align.PairwiseAlignmentResult{Int64,Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}}}:
  score: 11
  seq: 1 CGGATTA- 7
          || |||
  ref: 0 -GGTTTAC 7

```

`pairalign` takes an alignment type as its first argument, then two sequences to
align, and score (or cost) model. Alignment type is one of the following four
types:

* `GlobalAlignment`: global-to-global alignment
* `SemiGlobalAlignment`: local-to-global alignment
* `LocalAlignment`: local-to-local alignment
* `OverlapAlignment`: end-free alignment

For scoring model, `AffineGapScoreModel` is currently supported. It imposes an
**affine gap penalty** for insertions and deletions: `gap_open + k * gap_extend`
for a consecutive insertion/deletion of length `k`. The affine gap penalty is
flexible enough to create a constant and linear scoring model. Setting
`gap_extend = 0` or `gap_open = 0`, they are equivalent to the constant or
linear gap penalty, respectively.
The first argument of `AffineGapScoreModel` can be a substitution matrix like
`AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)`. For details on
substitution matrices, see the [Substitution matrices](@ref) section.

Alignment type can also be a distance of two sequences:

* `EditDistance`
* `LevenshteinDistance`
* `HammingDistance`

In this alignment, `CostModel` is used instead of `AffineGapScoreModel` to define
cost of substitution, insertion, and deletion:
```jlcon
julia> costmodel = CostModel(match=0, mismatch=1, insertion=1, deletion=1);

julia> pairalign(EditDistance(), "abcd", "adcde", costmodel)
Bio.Align.PairwiseAlignmentResult{Int64,String,String}:
  distance: 2
  seq: 1 abcd- 4
         | ||
  ref: 1 adcde 5

```

### Operations on pairwise alignment

`pairalign` returns a `PairwiseAlignmentResult` object and some accessors are
provided for it.

```@docs
score
distance
hasalignment
alignment
```

Pairwise alignment also implements some useful operations on it.

```@docs
count_matches
count_mismatches
count_insertions
count_deletions
count_aligned
```

The example below shows a use case of these operations:
```jlcon
julia> s1 = dna"CCTAGGAGGG";

julia> s2 = dna"ACCTGGTATGATAGCG";

julia> scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);

julia> res = pairalign(GlobalAlignment(), s1, s2, scoremodel)  # run pairwise alignment
Bio.Align.PairwiseAlignmentResult{Int64,Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}}}:
  score: 13
  seq:  0 -CCTAGG------AGGG 10
           ||| ||      || |
  ref:  1 ACCT-GGTATGATAGCG 16


julia> score(res)  # get the achieved score of this alignment
13

julia> aln = alignment(res)
PairwiseAlignment{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}}}:
  seq:  0 -CCTAGG------AGGG 10
           ||| ||      || |
  ref:  1 ACCT-GGTATGATAGCG 16


julia> count_matches(aln)
8

julia> count_mismatches(aln)
1

julia> count_insertions(aln)
1

julia> count_deletions(aln)
7

julia> count_aligned(aln)
17

julia> collect(aln)  # pairwise alignment is iterable
17-element Array{Tuple{Bio.Seq.DNA,Bio.Seq.DNA},1}:
 (-,A)
 (C,C)
 (C,C)
 (T,T)
 (A,-)
 (G,G)
 (G,G)
 (-,T)
 (-,A)
 (-,T)
 (-,G)
 (-,A)
 (-,T)
 (A,A)
 (G,G)
 (G,C)
 (G,G)

julia> DNASequence([x for (x, _) in aln])  # create aligned `s1` with gaps
17nt DNA Sequence:
-CCTAGG------AGGG

julia> DNASequence([y for (_, y) in aln])  # create aligned `s2` with gaps
17nt DNA Sequence:
ACCT-GGTATGATAGCG

```

## Substitution matrices

A substitution matrix is a function of substitution score (or cost) from one
symbol to other. Substitution value of `submat` from `x` to `y` can be obtained
by writing `submat[x,y]`.
In `Bio.Align`, `SubstitutionMatrix` and `DichotomousSubstitutionMatrix` are two
distinct types representing substitution matrices.

`SubstitutionMatrix` is a general substitution matrix type that is a thin
wrapper of regular matrix.

Some common substitution matrices are provided. For DNA and RNA, `EDNAFULL` is
defined:

```jlcon
julia> EDNAFULL
Bio.Align.SubstitutionMatrix{Bio.Seq.DNA,Int64}:
     A  C  G  T  M  R  W  S  Y  K  V  H  D  B  N
  A  5 -4 -4 -4  1  1  1 -4 -4 -4 -1 -1 -1 -4 -2
  C -4  5 -4 -4  1 -4 -4  1  1 -4 -1 -1 -4 -1 -2
  G -4 -4  5 -4 -4  1 -4  1 -4  1 -1 -4 -1 -1 -2
  T -4 -4 -4  5 -4 -4  1 -4  1  1 -4 -1 -1 -1 -2
  M  1  1 -4 -4 -1 -2 -2 -2 -2 -4 -1 -1 -3 -3 -1
  R  1 -4  1 -4 -2 -1 -2 -2 -4 -2 -1 -3 -1 -3 -1
  W  1 -4 -4  1 -2 -2 -1 -4 -2 -2 -3 -1 -1 -3 -1
  S -4  1  1 -4 -2 -2 -4 -1 -2 -2 -1 -3 -3 -1 -1
  Y -4  1 -4  1 -2 -4 -2 -2 -1 -2 -3 -1 -3 -1 -1
  K -4 -4  1  1 -4 -2 -2 -2 -2 -1 -3 -3 -1 -1 -1
  V -1 -1 -1 -4 -1 -1 -3 -1 -3 -3 -1 -2 -2 -2 -1
  H -1 -1 -4 -1 -1 -3 -1 -3 -1 -3 -2 -1 -2 -2 -1
  D -1 -4 -1 -1 -3 -1 -1 -3 -3 -1 -2 -2 -1 -2 -1
  B -4 -1 -1 -1 -3 -3 -3 -1 -1 -1 -2 -2 -2 -1 -1
  N -2 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
(underlined values are default ones)

```

For amino acids, PAM (Point Accepted Mutation) and BLOSUM (BLOcks SUbstitution Matrix) matrices are defined:

```jlcon
julia> BLOSUM62
Bio.Align.SubstitutionMatrix{Bio.Seq.AminoAcid,Int64}:
     A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  O  U  B  J  Z  X  *
  A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0  0̲  0̲ -2  0̲ -1  0 -4
  R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3  0̲  0̲ -1  0̲  0 -1 -4
  N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  0̲  0̲  3  0̲  0 -1 -4
  D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  0̲  0̲  4  0̲  1 -1 -4
  C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1  0̲  0̲ -3  0̲ -3 -2 -4
  Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0̲  0̲  0  0̲  3 -1 -4
  E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  0̲  0̲  1  0̲  4 -1 -4
  G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3  0̲  0̲ -1  0̲ -2 -1 -4
  H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0̲  0̲  0  0̲  0 -1 -4
  I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3  0̲  0̲ -3  0̲ -3 -1 -4
  L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1  0̲  0̲ -4  0̲ -3 -1 -4
  K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0̲  0̲  0  0̲  1 -1 -4
  M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1  0̲  0̲ -3  0̲ -1 -1 -4
  F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1  0̲  0̲ -3  0̲ -3 -1 -4
  P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2  0̲  0̲ -2  0̲ -1 -2 -4
  S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0̲  0̲  0  0̲  0  0 -4
  T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0  0̲  0̲ -1  0̲ -1  0 -4
  W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3  0̲  0̲ -4  0̲ -3 -2 -4
  Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1  0̲  0̲ -3  0̲ -2 -1 -4
  V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4  0̲  0̲ -3  0̲ -2 -1 -4
  O  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲
  U  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲
  B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  0̲  0̲  4  0̲  1 -1 -4
  J  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲  0̲
  Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  0̲  0̲  1  0̲  4 -1 -4
  X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1  0̲  0̲ -1  0̲ -1 -1 -4
  * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  0̲  0̲ -4  0̲ -4 -4  1
(underlined values are default ones)

```

| Matrix   | Constants                                                  |
| :------- | :----------                                                |
| PAM      | `PAM30`, `PAM70`, `PAM250`                                 |
| BLOSUM   | `BLOSUM45`, `BLOSUM50`, `BLOSUM62`, `BLOSUM80`, `BLOSUM90` |

These matrices are downloaded from: <ftp://ftp.ncbi.nih.gov/blast/matrices/>.

`SubstitutionMatrix` can be modified like a regular matrix:

```jlcon
julia> mysubmat = copy(BLOSUM62);  # create a copy

julia> mysubmat[AA_A,AA_R]  # score of AA_A => AA_R substitution is -1
-1

julia> mysubmat[AA_A,AA_R] = -3  # set the score to -3
-3

julia> mysubmat[AA_A,AA_R]  # the score is modified
-3

```

Make sure to create a copy of the original matrix when you create a matrix from
a predefined matrix. In the above case, `BLOSUM62` is shared in the whole
program and modification on it will affect any result that uses `BLOSUM62`.

`DichotomousSubstitutionMatrix` is a specialized matrix for matching or
mismatching substitution.  This is a preferable choice when performance is more
important than flexibility because looking up score is faster than
`SubstitutionMatrix`.

```jlcon
julia> submat = DichotomousSubstitutionMatrix(1, -1)
Bio.Align.DichotomousSubstitutionMatrix{Int64}:
     match =  1
  mismatch = -1

julia> submat['A','A']  # match
1

julia> submat['A','B']  # mismatch
-1

```


## Alignment file formats for high-throughput sequencing

High-throughput sequencing (HTS) technologies generate a large amount of data in
the form of a large number of nucleotide sequencing reads. One of the most
common tasks in bioinformatics is to align these reads against known reference
genomes, chromosomes, or contigs. The `Bio.Align` module provides several data
formats commonly used for this kind of task.


### SAM and BAM file formats

SAM and BAM are the most popular file formats and have the same reading and
writing interface as all other formats in Bio.jl (see [Reading and writing
data](@ref) section). A typical code iterating over all records in a file looks
like below:
```julia
# import the SAM and BAM module
using Bio.Align

# open a BAM file
reader = open(BAM.Reader, "data.bam")

# iterate over BAM records
for record in reader
    # `record` is a BAM.Record object
    if BAM.ismapped(record)
        # print mapped position
        println(BAM.refname(record), ':', BAM.position(record))
    end
end

# close the BAM file
close(reader)
```

Accessor functions are defined in `SAM` and `BAM` modules.  Lists of these
functions to `SAM.Record` and `BAM.Record` are described in [SAM](@ref) and
[BAM](@ref) sections, respectively.

`SAM.Reader` and `BAM.Reader` implement the `header` function, which returns a
`SAM.Header` object. This is conceptually a sequence of `SAM.MetaInfo` objects
corresponding to header lines that start with '@' markers. To select
`SAM.MetaInfo` records with a specific tag, you can use the `find` function:
```jlcon
julia> reader = open(SAM.Reader, "data.sam");

julia> find(header(reader), "SQ")
7-element Array{Bio.Align.SAM.MetaInfo,1}:
 Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=Chr1 LN=30427671
 Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=Chr2 LN=19698289
 Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=Chr3 LN=23459830
 Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=Chr4 LN=18585056
 Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=Chr5 LN=26975502
 Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=chloroplast LN=154478
 Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=mitochondria LN=366924

```

A `SAM.MetaInfo` object can be created as follows:
```jlcon
julia> SAM.MetaInfo("SQ", ["SN" => "chr1", "LN" => 1234])
Bio.Align.SAM.MetaInfo:
    tag: SQ
  value: SN=chr1 LN=1234

julia> SAM.MetaInfo("CO", "comment")
Bio.Align.SAM.MetaInfo:
    tag: CO
  value: comment

```


### Performance tips

The size of a BAM file is often extremely large. The iterator interface
mentioned above allocates an object for each record and that may be a bottleneck
of reading data from a BAM file. In-place reading reuses a preallocated object
for every record and less memory allocation happens in reading:
```julia
reader = open(BAM.Reader, "data.bam")
record = BAM.Record()
while !eof(reader)
    read!(reader, record)
    # do something
end
```

Accessing optional fields will results in type instability in Julia, which has a
significant negative impact on performance. If the user knows the type of a
value in advance, specifying it as a type annotation will alleviate the problem:
```julia
for record in open(BAM.Reader, "data.bam")
    nm = record["NM"]::UInt8
    # do something
end
```
