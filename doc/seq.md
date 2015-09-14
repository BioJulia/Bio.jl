---
part: Seq
order: 1
title: Nucleotide and amino acid sequences
modules: [Bio.Seq]
...


```{.julia hide="true"}
using Bio.Seq
```


The `Seq` module provides representations and tools for manipulating nucleotide
and amino acid sequences. Sequences in Bio.jl are more strictly typed than in
many other libraries. They are special purpose types rather than simply strings.
Bio.jl currently provides three distinct sequence types: `DNASequence`,
`RNASequence`, and `AminoAcidSequence`. Each is defined over a separate
alphabet, represented by types `DNANucleotide`, `RNANucleotide`, and
`AminoAcid`, respectively. Symbols from multiple alphabets can't be intermixed
in one sequence type.

Though this strictness sacrifices some convenience, it also means you can always
rely on a `DNASequence` to store DNA and nothing but DNA, without having to
check, or deal with lowercase versus uppercase and so on. Strict separation of
sequence types also means we are free to choose the most efficient
representation. DNA and RNA sequences are encoded using two bits per base making
them extremely memory efficient, and also allowing us to speed up many common
operations like nucleotide composition, reverse complement, and k-mer
enumeration.

Sequence are all able to represent missing or unobserved values using the
special symbols `DNA_N`, `RNA_N`, and `AA_X`.


# Constructing sequences and nucleotides

Nucleotide or amino acid symbols can be constructed by converting regular
characters.

```{.julia execute="false"}
# These calls produce three different types each representing "A" in a different alphabet
convert(AminoAcid, 'A')
convert(DNANucleotide, 'A')
convert(RNANucleotide, 'A')
```

There are also constants defined for each, such as `DNA_A`, `RNA_A`, and `AA_A`.

Sequence types corresponding to these alphabets can be constructed a number of
different ways. Most immediately, sequence literals can be constructed using
the string macros `dna`, `rna`, and `aa`.

```{.julia execute="false"}
# String decorators are provided for common sequence types
dna"TACGTANNATC"
rna"AUUUGNCCANU"
aa"ARNDCQEGHILKMFPSTWYVX"
```

Sequence can be constructed from strings or arrays of nucleotide or amino acid
symbols using constructors or the `convert` function.

```{.julia execute="false"}
DNASequence("TTANGTAGACCG")
DNASequence([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])
```

Using `convert`, these operations are reversible: sequences can be converted to
strings or arrays:

```{.julia execute="false"}
convert(String, dna"TTANGTAGACCG")
#convert(Vector{DNANucleotide}, dna"TTANGTAGACCG")
```

Sequences can also be concatenated into longer sequences

```{.julia execute="false"}
DNASequence(dna"ACGT", dna"NNNN", dna"TGCA")
```

Despite being separate types, `DNASequence` and `RNASequence` can freely be
converted between efficiently without copying the underlying data.

```{.julia execute="false"}
convert(RNASequence, dna"TTANGTAGACCG")
```

A translatable `RNASequence` can also be converted to an `AminoAcidSequence`
using the `translate` function described below.

# Indexing and iteration

Sequences for the most part behave like other string types. They can be indexed
using integers or ranges:

```{.julia execute="false"}
seq = dna"ACGTTTANAGTNNAGTACC"
seq[5]
seq[6:end]
```

They also work as iterators over nucleotides or amino acids:

```{.julia execute="false"}
for (i, nt) in enumerate(seq)
    if nt == DNA_N
        println("N in position ", i)
    end
end
```

# Mutability

Most sequences in Bio.jl are immutable, which means attempting to change
nucleotides in the sequence will result in an error. The advantage of
immutability is that subsequences don't have to make copies of the original data
but can instead just point to different intervals in the same memory.

```{.julia execute="false"}
# initialize a random million nucleotide sequence
large_sequence = DNASequence(rand([DNA_A, DNA_C, DNA_T, DNA_G], 1000000))

# This subsequence requires very little additional memory, since it doesn't copy the underlying data.
suffix = large_sequence[2:end]
```

Immutable sequences have an obvious disadvantage: some algorithms are easier or
more efficient to implement if a sequence can be altered or overwritten. For
this reason we provide a mutable mode for sequences. Mutability can be turned on
and off at any time by calling `mutable!(seq)` and `immutable!(seq)`,
respectively.

```{.julia execute="false"}
# This is an error, because sequences are immutable by default
seq = dna"ACGNACCTAGATAC"
seq[1] = DNA_T

# Switching the sequence into mutable mode allows this
mutable!(seq)
seq[1] = DNA_T
```

There are few subtleties to mode switching to be aware of if you're writing
performance critical code. Switching a mutable sequence to immutable is always
efficient; no data needs to be copied to do so. So is switching from immutable
to mutable, but only if the sequence hasn't had any subsequences made, and isn't
itself a subsequence, otherwise a complete copy of the data will be made.

# Operations on sequences

A number of common sequence operations on nucleotide sequences are provided in the Seq module.

{{reverse_complement}}


{{reverse}}


{{complement}}


{{repeat}}


{{mismatches}}


{{NucleotideCounts}}


{{translate}}

The `Seq` module contains all NCBI defined genetic codes:
```{.julia execute="false"}
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

In most cases, and by default, `standard_genetic_code` is used.

# Nucleotide K-mers

A common strategy to simplify the analysis of sequence data is to operate or
short k-mers, for size fixed size `k`. These can be packed into machine integers
allowing extremely efficient code. The Seq module has built in support for
representing short sequences in 64-bit integers. Besides being fixed length,
`Kmer` types, unlike other sequence types cannot contain `N` symbols.

The `Kmer{T, K}` type parameterized on alphabet (`T`, either `DNANucleotide`, or
`RNANucleotide`) and size `K`. A number of functions are provided for operating
on `Kmers`.

{{each}}

{{KmerCounts}}

{{canonical}}

{{neighbors}}


# Sequence records

The `SeqRecord` type is used to represent a named sequence, optionally with
accompanying metadata. It is defined as:
sequence. It is as:
```{.julia execute="false"}
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

