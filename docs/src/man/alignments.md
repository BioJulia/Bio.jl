# Align: Sequence Alignments

    {meta}
    CurrentModule = Bio.Align
    DocTestSetup = quote
        using Bio.Seq
        using Bio.Align
    end

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
```
              4   8   12   17    20 23  27
              |   |   |    |     |  |   |
reference: AGGGTGGCATTTATCAG---ACGTTTCGAGAC
    query:     TGGC----ATCATTTAACG---CAAG
              |   |        |  |  |      |
              0   4        9  12 15     19

```

Using anchors we would represent this as the following series of anchors:
```
AlignmentAnchor(0, 4, OP_START),
AlignmentAnchor(4, 8, OP_MATCH),
AlignmentAnchor(4, 12, OP_DELETE),
AlignmentAnchor(9, 17, OP_MATCH),
AlignmentAnchor(12, 17, OP_INSERT),
AlignmentAnchor(15, 20, OP_MATCH),
AlignmentAnchor(15, 23, OP_DELETE),
AlignmentAnchor(19, 27, OP_MATCH)
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


## Operating on Alignments

    {docs}
    first
    last
    cigar


## Pairwise Alignment

Pairwise alignment is a sequence alignment between two sequences.
The `Bio.Align` module provides some pairwise alignment algorithms that
optimize alignment score.

Let $\Sigma$ be a set of symbols, and $a$ and $b$ be two sequences of $\Sigma$, then
pairwise alignment is defined as a pair of two gap-padded sequences $a^\prime$
and $b^\prime$ of $\Sigma \cap \\{\texttt{'-'}\\}$. Alignment score of this
alignment is defined as a sum of element-wise scores:
$$ \sum_i score(a^\prime_i, b^\prime_i) $$
where $score(x, y)$ is a scoring function of two symbols $x$ and $y$.

A pair of optimization problem and score (or cost) model determines the best
alignment between two sequences. For example, a pair of `GlobalAlignment` problem and
`AffineGapScoreModel` model is well known as end-to-end alignment with affine
penalty for gap insertion. Two strings, "succeed" and "precede", are aligned as
follows in this settings:
```julia
julia> using Bio.Align

julia> problem = GlobalAlignment()
Bio.Align.GlobalAlignment()

julia> model = AffineGapScoreModel(
           match=5,
           mismatch=-2,
           gap_open=-3,
           gap_extend=-1
       )
Bio.Align.AffineGapScoreModel{Int64}:
       match = 5
    mismatch = -2
    gap_open = -3
  gap_extend = -1

julia> pairalign(problem, "succeed", "precede", model)
PairwiseAlignmentResult{Int64,ASCIIString,ASCIIString}:
  score: 1
  seq: succe-ed
          || |
  ref: precede-

```


### Alignment types

#### Alignments

| Alignment Problem | Description |
|:---------------|:-----------|
| `GlobalAlignment` | global-global alignment |
| `SemiGlobalAlignment` | global-local alignment |
| `OverlapAlignment` | global-global alignment without end gap penalties |
| `LocalAlignment` | local-local alignment |

#### Global vs. local alignment

```julia
julia> a = dna"CAGGTAGTAGAGTATATTATGGCCATTTCTATCGTTATTT"
40nt DNA Sequence:
CAGGTAGTAGAGTATATTATGGCCATTTCTATCGTTATTT

julia> b = dna"ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA"
40nt DNA Sequence:
ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA

julia> affinegap = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-10, gap_extend=-1)
Bio.Align.AffineGapScoreModel{Int64}:
       match = 5
    mismatch = -4
    gap_open = -10
  gap_extend = -1

julia> pairalign(GlobalAlignment(), a, b, affinegap)
Bio.Seq.PairwiseAlignmentResult{Int64,Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Bio.Seq.BioSeque
nce{Bio.Seq.DNAAlphabet{4}}}:
  score: 12
  seq:  0 ----CAGGTAGTAGAGTATATTATGGCCATT---TCTATCGTTATTT 40
              |   || |    ||||| |||| ||||   |||  |   |
  ref:  1 ATCTCTAATATT----TATATCATGGACATTAAATCTCGCAGGA--- 40


julia> pairalign(LocalAlignment(), a, b, affinegap)
Bio.Seq.PairwiseAlignmentResult{Int64,Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Bio.Seq.BioSeque
nce{Bio.Seq.DNAAlphabet{4}}}:
  score: 59
  seq: 13 TATATTATGGCCATT---TCT 30
          ||||| |||| ||||   |||
  ref: 13 TATATCATGGACATTAAATCT 33


```

#### Global vs. semi-global alignment
```
julia> a = dna"TATTCCATGATT"
12nt DNA Sequence:
TATTCCATGATT

julia> b = dna"ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA"
40nt DNA Sequence:
ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA

julia> affinegap = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-10, gap_extend=-1)
Bio.Align.AffineGapScoreModel{Int64}:
       match = 5
    mismatch = -4
    gap_open = -10
  gap_extend = -1

julia> pairalign(GlobalAlignment(), a, b, affinegap)
Bio.Align.PairwiseAlignmentResult{Int64,Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Bio.Se
q.BioSequence{Bio.Seq.DNAAlphabet{4}}}:
  score: -16
  seq:  0 ------------TATTCCATG---ATT------------- 12
                      |||  ||||   |||
  ref:  1 ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA 40


julia> pairalign(SemiGlobalAlignment(), a, b, affinegap)
Bio.Align.PairwiseAlignmentResult{Int64,Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Bio.Se
q.BioSequence{Bio.Seq.DNAAlphabet{4}}}:
  score: 29
  seq:  0 ------------TATTCCATG---ATT------------- 12
                      |||  ||||   |||
  ref:  1 ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA 40


```

#### Global vs. overlap alignment

```
julia> a = dna"ATCTAACATTGGACATTAAATCTCGCATGATCGGACATTG"

julia> b = dna"ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA"
40nt DNA Sequence:
ATCTCTAATATTTATATCATGGACATTAAATCTCGCAGGA

julia> affinegap = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-10, gap_extend=-1)
Bio.Align.AffineGapScoreModel{Int64}:
       match = 5
    mismatch = -4
    gap_open = -10
  gap_extend = -1

```

### Distances

| Alignment Problem | Description |
|:----------------|:-------------|
| `EditDistance` | edit distance |
| `LevenshteinDistance` | special case of edit distance |
| `HammingDistance` | special case of edit distance |


### Scoring models


### Substitution matrices

A substitution matrix is a function of substitution score (or cost) from one
symbol to other. Substitution value of `submat` from `x` to `y` can be obtained
by writing `submat[x,y]`.
In `Bio.Align`, `SubstitutionMatrix` and `DichotomousSubstitutionMatrix` are two
distinct types representing substitution matrices.

`SubstitutionMatrix` is a general substitution matrix type that is a thin
wrapper of regular matrix.

Some common substitution matrices are provided. For DNA and RNA, `EDNAFULL` is
defined:

```julia
julia> EDNAFULL
Bio.Align.SubstitutionMatrix{Bio.Seq.DNANucleotide,Int64}:
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

```julia
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

```julia
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

```julia
julia> submat = DichotomousSubstitutionMatrix(1, -1)
Bio.Align.DichotomousSubstitutionMatrix{Int64}:
     match =  1
  mismatch = -1

julia> submat['A','A']  # match
1

julia> submat['A','B']  # mismatch
-1

```
