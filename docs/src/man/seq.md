# Bio.Seq: Biological Sequences

    {meta}
    CurrentModule = Bio.Seq
    DocTestSetup = quote
        using Bio.Seq
    end

The `Seq` module provides representations and tools for manipulating nucleotide
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
the string macros `dna`, `rna`, and `aa`:

```julia
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

```

Sequence can be constructed from strings or arrays of nucleotide or amino acid
symbols using constructors or the `convert` function:

```julia
julia> DNASequence("TTANC")
5nt DNA Sequence:
TTANC

julia> DNASequence([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
5nt DNA Sequence:
TTANC

```

Using `convert`, these operations are reversible: sequences can be converted to
strings or arrays:

```julia
julia> convert(ASCIIString, dna"TTANGTA")
"TTANGTA"

julia> convert(Vector{DNANucleotide}, dna"TTANGTA")
7-element Array{Bio.Seq.DNANucleotide,1}:
 T
 T
 A
 N
 G
 T
 A

```

Sequences can also be concatenated into longer sequences:

```julia
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

```julia
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
using the `translate` function described below.


## Indexing and modifying

Sequences for the most part behave like other vector or string types. They can be indexed
using integers or ranges:

```julia
julia> seq = dna"ACGTTTANAGTNNAGTACC"
19nt DNA Sequence:
ACGTTTANAGTNNAGTACC

julia> seq[5]
T

julia> seq[6:end]
14nt DNA Sequence:
TANAGTNNAGTACC

```

Indexing by range creates a subsequence of the original sequence. Unlike
`ASCIIString` and `Vector` in the standard library, creating a subsequences is
copy-free: a subsequence is just a reference to the original sequence with its range.
You may think that this is unsafe because modifying subsequences propagates to
the original sequence, but this doesn't happen actually:

```julia
julia> seq = dna"AAAA"  # create a sequence
4nt DNA Sequence:
AAAA

julia> subseq = seq[1:2]  # create a subsequence from `seq`
2nt DNA Sequence:
AA

julia> subseq[2] = DNA_T  # modify the second element of it
T

julia> subseq  # the subsequence is modified
2nt DNA Sequence:
AT

julia> seq  # but the original sequence is not
4nt DNA Sequence:
AAAA

```

This is because modifying a sequence checks whether its underlying data are shared with other sequences under the hood. If and only if the data are shared, the subsequence creates a copy of itself. Any modifying operation does this check. This is called *copy-on-write* strategy and users don't need to care about it because it is transparent from outward.

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

    {docs}
    complement!(seq)
    reverse_complement!(seq)


```julia
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

julia> complement!(seq)
6nt DNA Sequence:
ATACGT

julia> reverse_complement!(seq)
6nt DNA Sequence:
ACGTAT

```

Sequences also work as iterators over symbols:

```julia
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

    {docs}
    complement
    reverse_complement
    mismatches
    composition
    translate

The `Bio.Seq` module contains all NCBI defined genetic codes:

```julia
standard_genetic_code,
vertebrate_mitochondrial_genetic_code,
yeast_mitochondrial_genetic_code,
mold_mitochondrial_genetic_code,
invertebrate_mitochondrial_genetic_code,
ciliate_nuclear_genetic_code,
echinoderm_mitochondrial_genetic_code,
euplotid_nuclear_genetic_code,
bacterial_plastid_genetic_code,
alternative_yeast_nuclear_genetic_code,
ascidian_mitochondrial_genetic_code,
alternative_flatworm_mitochondrial_genetic_code,
chlorophycean_mitochondrial_genetic_code,
trematode_mitochondrial_genetic_code,
scenedesmus_obliquus_mitochondrial_genetic_code,
thraustochytrium_mitochondrial_genetic_code,
pterobrachia_mitochondrial_genetic_code,
candidate_division_sr1_genetic_code
```

<http://www.insdc.org/documents/feature_table.html#7.4.5>

In most cases, and by default, `standard_genetic_code` is used.


## Compact representation

As we saw above, DNA and RNA sequences can store any ambiguous nucleotides like 'N'.
If you are sure that nucleotide sequences store unambiguous nucleotides only, you can
save the memory space sequences use. `DNAAlphabet{2}` is an alphabet that uses
two bits per base and limited only to umambiguous nucleotide symbols, namely ACGT in DNA and ACGU in RNA. To create a sequence of this alphabet, you need to
explicitly pass `DNAAlphabet{2}` to `BioSequence` as its parametric type:

```julia
julia> seq = BioSequence{DNAAlphabet{2}}("ACGT")
4nt DNA Sequence:
ACGT

```

Recall that `DNASequence` is a type alias of `BioSequence{DNAAlphabet{4}}`,
which uses four bits per base. That is, `BioSequence{DNAAlphabet{2}}` saves
half memory footprint of `BioSequence{DNAAlphabet{4}}`. If you need to handle
reference sequences that are composed of five nucleotides, ACGTN, consider to use
[ReferenceSequences.jl](https://github.com/BioJulia/ReferenceSequences.jl). This
compresses positions of 'N' and enables to handle long DNA sequences with the
near space of two-bit encoding.


## Nucleotide K-mers

A common strategy to simplify the analysis of sequence data is to operate or
short k-mers, for size fixed size `k`. These can be packed into machine integers
allowing extremely efficient code. The `Bio.Seq` module has built in support for
representing short sequences in 64-bit integers. Besides being fixed length,
`Kmer` types, unlike other sequence types cannot contain ambiguous symbols like
'N'.

The `Kmer{T,k}` type parameterized on alphabet (`T`, either `DNANucleotide`, or
`RNANucleotide`) and size `k`. A number of functions are provided for operating
on `Kmers`.

    {docs}
    each
    KmerCounts
    canonical
    neighbors


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

