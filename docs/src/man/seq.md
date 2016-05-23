# Bio.Seq: Biological Sequences

```@meta
CurrentModule = Bio.Seq
DocTestSetup = quote
    using Bio.Seq
end
```

The `Bio.Seq` module provides several data types for handling biological symbols
and sequences.


## Biological symbols

```@meta
CurrentModule = Bio.Seq
DocTestSetup = quote
    using Bio.Seq
end
```

The `Bio.Seq` module provides three biological symbol (character) types:

| Type            | Meaning        |
| :-------------- | :------------- |
| `DNANucleotide` | DNA nucleotide |
| `RNANucleotide` | RNA nucleotide |
| `AminoAcid`     | Amino acid     |

These symbols can be elements of sequences like characters can be elements of
strings.  See [Sequences](@ref) and [Nucleotide k-mers](@ref) section for
details.


### DNA and RNA nucleotides

Set of nucleotide symbols in Bio.jl covers IUPAC nucleotide base plus a gap symbol:

| Symbol | Constant              | Meaning                  |
| :----- | :-------------------- | :----------------------- |
| 'A'    | `DNA_A` / `RNA_A`     | A; Adenine               |
| 'C'    | `DNA_C` / `RNA_C`     | C; Cytosine              |
| 'G'    | `DNA_G` / `RNA_G`     | G; Guanine               |
| 'T'    | `DNA_T`               | T; Thymine (DNA only)    |
| 'U'    | `RNA_U`               | U; Uracil (RNA only)     |
| 'M'    | `DNA_M` / `RNA_M`     | A or C                   |
| 'R'    | `DNA_R` / `RNA_R`     | A or G                   |
| 'W'    | `DNA_W` / `RNA_W`     | A or T                   |
| 'S'    | `DNA_S` / `RNA_S`     | C or G                   |
| 'Y'    | `DNA_Y` / `RNA_Y`     | C or T                   |
| 'K'    | `DNA_K` / `RNA_K`     | G or T                   |
| 'V'    | `DNA_V` / `RNA_V`     | A or C or G; not T       |
| 'H'    | `DNA_H` / `RNA_H`     | A or C or T; not G       |
| 'D'    | `DNA_D` / `RNA_D`     | A or G or T; not C       |
| 'B'    | `DNA_B` / `RNA_B`     | C or G or T; not A       |
| 'N'    | `DNA_N` / `RNA_N`     | A or C or G or T         |
| '-'    | `DNA_Gap` / `RNA_Gap` | Gap (none of the above)  |

<http://www.insdc.org/documents/feature_table.html#7.4.1>

Symbols are accessible as constants with `DNA_` or `RNA_` prefix:

```jlcon
julia> DNA_A
DNA_A

julia> DNA_T
DNA_T

julia> RNA_U
RNA_U

julia> DNA_Gap
DNA_Gap

julia> typeof(DNA_A)
Bio.Seq.DNANucleotide

julia> typeof(RNA_A)
Bio.Seq.RNANucleotide

```

Symbols can be constructed by converting regular characters:

```jlcon
julia> convert(DNANucleotide, 'C')
DNA_C

julia> convert(DNANucleotide, 'C') === DNA_C
true

```


### Amino acids

Set of amino acid symbols also covers IUPAC amino acid symbols plus a gap symbol:

| Symbol       | Constant        | Meaning                     |
| :----------- | :-------------- | :-------------------------- |
| 'A'          | `AA_A`          | Alanine                     |
| 'R'          | `AA_R`          | Arginine                    |
| 'N'          | `AA_N`          | Asparagine                  |
| 'D'          | `AA_D`          | Aspartic acid (Aspartate)   |
| 'C'          | `AA_C`          | Cysteine                    |
| 'Q'          | `AA_Q`          | Glutamine                   |
| 'E'          | `AA_E`          | Glutamic acid (Glutamate)   |
| 'G'          | `AA_G`          | Glycine                     |
| 'H'          | `AA_H`          | Histidine                   |
| 'I'          | `AA_I`          | Isoleucine                  |
| 'L'          | `AA_L`          | Leucine                     |
| 'K'          | `AA_K`          | Lysine                      |
| 'M'          | `AA_M`          | Methionine                  |
| 'F'          | `AA_F`          | Phenylalanine               |
| 'P'          | `AA_P`          | Proline                     |
| 'S'          | `AA_S`          | Serine                      |
| 'T'          | `AA_T`          | Threonine                   |
| 'W'          | `AA_W`          | Tryptophan                  |
| 'Y'          | `AA_Y`          | Tyrosine                    |
| 'V'          | `AA_V`          | Valine                      |
| 'O'          | `AA_O`          | Pyrrolysine                 |
| 'U'          | `AA_U`          | Selenocysteine              |
| 'B'          | `AA_B`          | Aspartic acid or Asparagine |
| 'J'          | `AA_J`          | Leucine or Isoleucine       |
| 'Z'          | `AA_Z`          | Glutamine or Glutamic acid  |
| 'X'          | `AA_X`          | Any amino acid              |
| '*'          | `AA_Term`       | Termination codon           |
| '-'          | `AA_Gap`        | Gap (none of the above)     |

<http://www.insdc.org/documents/feature_table.html#7.4.3>

Symbols are accessible as constants with `AA_` prefix:
```jlcon
julia> AA_A
AA_A

julia> AA_Q
AA_Q

julia> AA_Term
AA_Term

julia> typeof(AA_A)
Bio.Seq.AminoAcid

```

Symbols can be constructed by converting regular characters:
```jlcon
julia> convert(AminoAcid, 'A')
AA_A

julia> convert(AminoAcid, 'P') === AA_P
true

```


### Arithmetic

Biological symbols behaves like `Char`:
```jlcon
julia> DNA_A == DNA_A  # equivalence
true

julia> DNA_A < DNA_C < DNA_G < DNA_T  # order
true

julia> DNA_A + 1  # addition
DNA_C

julia> DNA_T - 3  # subtraction
DNA_A

julia> DNA_T - DNA_C  # difference
2

```

Note that these operations are cyclic:
```jlcon
julia> DNA_C + 15
DNA_A

julia> DNA_A - 1
DNA_Gap

```


### Symbol Ranges

Consecutive symbol sets can be created using a colon like integer ranges:
```jlcon
julia> DNA_A:DNA_T    # unambiguous DNA nucleotides (A, C, G, T)
DNA_A:DNA_T

julia> DNA_A:DNA_N    # all DNA nucleotides except gap
DNA_A:DNA_N

julia> DNA_A:DNA_Gap  # all DNA nucleotides including gap
DNA_A:DNA_Gap

julia> AA_A:AA_V      # standard amino acids
AA_A:AA_V

julia> AA_A:AA_U      # standard amino acids + pyrrolysine (O) + selenocysteine (U)
AA_A:AA_U

julia> AA_A:AA_X      # all amino acids except teminal codon (*) and gap
AA_A:AA_X

julia> AA_A:AA_Gap    # all amino acids including teminal codon (*) and gap
AA_A:AA_Gap

```

Most range operations are supported, especially the iterator interface and the
membership operator (`in` or `∈`) will be useful in many situations:
```jlcon
julia> for nt in DNA_A:DNA_T; println(nt); end
A
C
G
T

julia> DNA_C in DNA_A:DNA_T  # DNA_C is in the range of unambiguous DNA nucleotides
true

julia> DNA_N in DNA_A:DNA_T  # DNA_N is not in it
false

```


### Other functions

```@docs
alphabet
gap
iscompatible
```


## Sequences

The `Bio.Seq` module provides representations and tools for manipulating nucleotide
and amino acid sequences. Sequences in Bio.jl are more strictly typed than in
many other libraries. They are special purpose types rather than simply strings.
Bio.jl currently provides a single generic sequence type:
`BioSequence{A<:Alphabet}`. `BioSequence` is parameterized by an alphabet type
`A` that defines the domain (or set) of biological symbols, and each alphabet
has an associated symbol type. For example, `AminoAcidAlphabet` is associated
with `AminoAcid` and hence an object of `BioSequence{AminoAcidAlphabet}`
represents a sequence of amino acids.  Symbols from multiple alphabets can't
be intermixed in one sequence type.

The following table summarizes common sequence types that are defined in the
`Bio.Seq` module:

| Type                               | Symbol type      | Type alias          |
| :--------------------------------- | :--------------- | :------------------ |
| `BioSequence{DNAAlphabet{4}}`      | `DNANucleotide`  | `DNASequence`       |
| `BioSequence{RNAAlphabet{4}}`      | `RNANucleotide`  | `RNASequence`       |
| `BioSequence{AminoAcidAlphabet}`   | `AminoAcid`      | `AminoAcidSequence` |
| `BioSequence{CharAlphabet}`        | `Char`           | `CharSequence`      |


Though this strictness sacrifices some convenience, it also means you can always
rely on a `DNASequence` to store DNA and nothing but DNA, without having to
check, or deal with lowercase versus uppercase and so on. Strict separation of
sequence types also means we are free to choose the most efficient
representation. DNA and RNA sequences are encoded using four bits per base by
default making them memory efficient, and also allowing us to speed up many
common operations like nucleotide composition, reverse complement, and *k*-mer
enumeration.


## Constructing sequences

Sequence types corresponding to these alphabets can be constructed a number of
different ways. Most immediately, sequence literals can be constructed using
the string macros `dna`, `rna`, `aa`, and `char`:
```jlcon
# String decorators are provided for common sequence types
julia> dna"TACGTANNATC"
11nt DNA Sequence:
TACGTANNATC

julia> rna"AUUUGNCCANU"
11nt RNA Sequence:
AUUUGNCCANU

julia> aa"ARNDCQEGHILKMFPSTWYVX"
21aa Amino Acid Sequence:
ARNDCQEGHILKMFPSTWYVX

julia> char"αβγδϵ"
5char Char Sequence:
αβγδϵ

```

Sequence can also be constructed from strings or arrays of nucleotide or amino
acid symbols using constructors or the `convert` function:
```jlcon
julia> DNASequence("TTANC")
5nt DNA Sequence:
TTANC

julia> DNASequence([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
5nt DNA Sequence:
TTANC

julia> convert(DNASequence, [DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
5nt DNA Sequence:
TTANC

```

Using `convert`, these operations are reversible: sequences can be converted to
strings or arrays:
```jlcon
julia> convert(ASCIIString, dna"TTANGTA")
"TTANGTA"

julia> convert(Vector{DNANucleotide}, dna"TTANGTA")
7-element Array{Bio.Seq.DNANucleotide,1}:
 DNA_T
 DNA_T
 DNA_A
 DNA_N
 DNA_G
 DNA_T
 DNA_A

```

Sequences can also be concatenated into longer sequences:
```jlcon
julia> DNASequence(dna"ACGT", dna"NNNN", dna"TGCA")
12nt DNA Sequence:
ACGTNNNNTGCA

julia> dna"ACGT" * dna"TGCA"
8nt DNA Sequence:
ACGTTGCA

julia> repeat(dna"TA", 10)
20nt DNA Sequence:
TATATATATATATATATATA

julia> dna"TA" ^ 10
20nt DNA Sequence:
TATATATATATATATATATA

```

Despite being separate types, `DNASequence` and `RNASequence` can freely be
converted between efficiently without copying the underlying data:
```jlcon
julia> dna = dna"TTANGTAGACCG"
12nt DNA Sequence:
TTANGTAGACCG

julia> rna = convert(RNASequence, dna)
12nt RNA Sequence:
UUANGUAGACCG

julia> dna.data === rna.data  # underlying data are same
true

```

A translatable `RNASequence` can also be converted to an `AminoAcidSequence`
using the [`translate`](@ref) function described below.


## Indexing and modifying

Sequences for the most part behave like other vector or string types. They can
be indexed using integers or ranges:
```jlcon
julia> seq = dna"ACGTTTANAGTNNAGTACC"
19nt DNA Sequence:
ACGTTTANAGTNNAGTACC

julia> seq[5]
DNA_T

julia> seq[6:end]
14nt DNA Sequence:
TANAGTNNAGTACC

```

Indexing by range creates a subsequence of the original sequence. Unlike
`ASCIIString` and `Vector` in the standard library, creating a subsequences is
copy-free: a subsequence is just a reference to the original sequence with its
range.  You may think that this is unsafe because modifying subsequences
propagates to the original sequence, but this doesn't happen actually:
```jlcon
julia> seq = dna"AAAA"  # create a sequence
4nt DNA Sequence:
AAAA

julia> subseq = seq[1:2]  # create a subsequence from `seq`
2nt DNA Sequence:
AA

julia> subseq[2] = DNA_T  # modify the second element of it
DNA_T

julia> subseq  # the subsequence is modified
2nt DNA Sequence:
AT

julia> seq  # but the original sequence is not
4nt DNA Sequence:
AAAA

```

This is because modifying a sequence checks whether its underlying data are
shared with other sequences under the hood. If and only if the data are shared,
the subsequence creates a copy of itself. Any modifying operation does this
check. This is called *copy-on-write* strategy and users don't need to care
about it because it is transparent from outward.

The following modifying operations are currently supported:
```julia
setindex!(seq, item, index)
push!(seq, item)
pop!(seq)
shift!(seq)
unshift!(seq, item)
insert!(seq, index, item)
deleteat!(seq, index)
append!(seq, other_seq)
copy!(dst_seq, dest_offset, src_seq, src_offset, len)
reverse!(seq)
```

```@docs
complement!(seq)
reverse_complement!(seq)
```


```jlcon
julia> seq = dna"ACG"
3nt DNA Sequence:
ACG

julia> push!(seq, DNA_T)
4nt DNA Sequence:
ACGT

julia> append!(seq, dna"AT")
6nt DNA Sequence:
ACGTAT

julia> reverse!(seq)
6nt DNA Sequence:
TATGCA

julia> Seq.complement!(seq)  # semantically differs from Base.complement!
6nt DNA Sequence:
ATACGT

julia> reverse_complement!(seq)
6nt DNA Sequence:
ACGTAT

```

Sequences also work as iterators over symbols:
```jlcon
julia> n = 0
0

julia> for nt in dna"ATNGNNT"
           if nt == DNA_N
               n += 1
           end
       end

julia> n
3

```


## Other operations on sequences

A number of common sequence operations are provided in the `Bio.Seq` module:

```@docs
complement
reverse_complement
mismatches
composition
```


### Translation

The [`translate`](@ref) funtion translates a sequence of codons in a RNA sequence
to a amino acid sequence besed on a genetic code mapping.  The `Bio.Seq` module
contains all NCBI defined genetic codes and they are registered in
[`ncbi_trans_table`](@ref).

```@docs
translate
ncbi_trans_table
```

```julia
julia> ncbi_trans_table
Translation Tables:
  1. The Standard Code (standard_genetic_code)
  2. The Vertebrate Mitochondrial Code (vertebrate_mitochondrial_genetic_code)
  3. The Yeast Mitochondrial Code (yeast_mitochondrial_genetic_code)
  4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (mold_mitochondrial_genetic_code)
  5. The Invertebrate Mitochondrial Code (invertebrate_mitochondrial_genetic_code)
  6. The Ciliate, Dasycladacean and Hexamita Nuclear Code (ciliate_nuclear_genetic_code)
  9. The Echinoderm and Flatworm Mitochondrial Code (echinoderm_mitochondrial_genetic_code)
 10. The Euplotid Nuclear Code (euplotid_nuclear_genetic_code)
 11. The Bacterial, Archaeal and Plant Plastid Code (bacterial_plastid_genetic_code)
 12. The Alternative Yeast Nuclear Code (alternative_yeast_nuclear_genetic_code)
 13. The Ascidian Mitochondrial Code (ascidian_mitochondrial_genetic_code)
 14. The Alternative Flatworm Mitochondrial Code (alternative_flatworm_mitochondrial_genetic_code)
 16. Chlorophycean Mitochondrial Code (chlorophycean_mitochondrial_genetic_code)
 21. Trematode Mitochondrial Code (trematode_mitochondrial_genetic_code)
 22. Scenedesmus obliquus Mitochondrial Code (scenedesmus_obliquus_mitochondrial_genetic_code)
 23. Thraustochytrium Mitochondrial Code (thraustochytrium_mitochondrial_genetic_code)
 24. Pterobranchia Mitochondrial Code (pterobrachia_mitochondrial_genetic_code)
 25. Candidate Division SR1 and Gracilibacteria Code (candidate_division_sr1_genetic_code)

```

<http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes>


## Compact representation

As we saw above, DNA and RNA sequences can store any ambiguous nucleotides like
'N'.  If you are sure that nucleotide sequences store unambiguous nucleotides
only, you can save the memory space of sequences. `DNAAlphabet{2}` is an
alphabet that uses two bits per base and limits to only unambiguous nucleotide
symbols (ACGT in DNA and ACGU in RNA). To create a sequence of this
alphabet, you need to explicitly pass `DNAAlphabet{2}` to `BioSequence` as its
parametric type:
```jlcon
julia> seq = BioSequence{DNAAlphabet{2}}("ACGT")
4nt DNA Sequence:
ACGT

```

Recall that `DNASequence` is a type alias of `BioSequence{DNAAlphabet{4}}`,
which uses four bits per base. That is, `BioSequence{DNAAlphabet{2}}` saves half
of memory footprint compared to `BioSequence{DNAAlphabet{4}}`. If you need to
handle reference sequences that are composed of five nucleotides, ACGTN,
consider to use
[ReferenceSequences.jl](https://github.com/BioJulia/ReferenceSequences.jl).
`ReferenceSequence`, exported from this package, compresses positions of 'N' and
enables to handle long DNA sequences with the near space of two-bit encoding.


## Nucleotide k-mers

A common strategy to simplify the analysis of sequence data is to operate or
short k-mers, for size fixed size `k`. These can be packed into machine integers
allowing extremely efficient code. The `Bio.Seq` module has built in support for
representing short sequences in 64-bit integers. Besides being fixed length,
`Kmer` types, unlike other sequence types cannot contain ambiguous symbols like
'N'.

The `Kmer{T,k}` type parameterized on alphabet (`T`, either `DNANucleotide`, or
`RNANucleotide`) and size `k`. A number of functions are provided for operating
on `Kmers`.

```@docs
each
KmerCounts
canonical
neighbors
```


## Sequence search

Several kinds of on-line search functions are provided.

### Exact search

Exact search functions search for an occurrence of the query symbol or
sequence. Four functions, `search`, `searchindex`, `rsearch`, and
`rsearchindex` are available:
```jlcon
julia> seq = dna"ACAGCGTAGCT";

julia> search(seq, DNA_G)  # search a query symbol
4:4

julia> query = dna"AGC";

julia> search(seq, query)  # search a query sequence
3:5

julia> searchindex(seq, query)
3

julia> rsearch(seq, query)  # similar to `search` but in the reverse direction
8:10

julia> rsearchindex(seq, query)  # similar to `searchindex` but in the reverse direction
8

```

These search functions take ambiguous symbols into account. That is, if two
symbols are compatible (e.g. `DNA_A` and `DNA_N`), they match when searching an
occurrence. In the following example, 'N' is a wild card that matches any
symbols:
```jlcon
julia> search(dna"ACNT", DNA_N)  # 'A' matches 'N'
1:1

julia> search(dna"ACNT", dna"CGT")  # 'N' matches 'G'
2:4

julia> search(dna"ACGT", dna"CNT")  # 'G' matches 'N'
2:4

```

The exact sequence search needs preprocessing phase of query sequence before
searching phase. This would be enough fast for most search applications. But
when searching a query sequence to large amounts of target sequences, caching
the result of preprocessing may save time. The `ExactSearchQuery` creates such
a preprocessed query object and is applicable to the search functions:
```jlcon
julia> query = ExactSearchQuery(dna"ATT");

julia> search(dna"ATTTATT", query)
1:3

julia> rsearch(dna"ATTTATT", query)
5:7

```


### Approximate search

The approximate search is similar to the exact search but allows a specific
number of errors. That is, it tries to find a subsequence of the target sequence
within a specific [Levenshtein
distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of the query
sequence:
```jlcon
julia> seq = dna"ACAGCGTAGCT";

julia> approxsearch(seq, dna"AGGG", 0)  # nothing matches with no errors
0:-1

julia> approxsearch(seq, dna"AGGG", 1)  # seq[3:5] matches with one error
3:6

julia> approxsearch(seq, dna"AGGG", 2)  # seq[1:4] matches with two errors
1:4

```

Like the exact search functions, four kinds of functions (`approxsearch`,
`approxsearchindex`, `approxrsearch`, and `approxrsearchindex`) are available:
```jlcon
julia> seq = dna"ACAGCGTAGCT"; pat = dna"AGGG";

julia> approxsearch(seq, pat, 2)        # return the range (forward)
1:4

julia> approxsearchindex(seq, pat, 2)   # return the starting index (forward)
1

julia> approxrsearch(seq, pat, 2)       # return the range (backward)
8:11

julia> approxrsearchindex(seq, pat, 2)  # return the starting index (backward)
8

```

Preprocessing can be cached in an `ApproximateSearchQuery` object:
```jlcon
julia> query = ApproximateSearchQuery(dna"AGGG");

julia> approxsearch(dna"AAGAGG", query, 1)
2:5

julia> approxsearch(dna"ACTACGT", query, 2)
4:6

```


## Sequence records

The `SeqRecord` type is used to represent a named sequence, optionally with
accompanying metadata. It is defined as:
sequence. It is as:
```julia
type SeqRecord{S, T}
    name::String
    seq::S
    metadata::T
end
```

The type of the `metadata` field depends on the source of the sequence record.
For example, if a record is read from a FASTA file, metadata contains the
description field. If from a FASTQ file, a quality scores assigned to base calls
during sequencing.
