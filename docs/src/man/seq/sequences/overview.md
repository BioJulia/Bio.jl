# Bioloigcal sequence data-types

The `Bio.Seq` module provides representations and tools for manipulating
nucleotide and amino acid sequences. Sequences in Bio.jl are more strictly typed
than in many other libraries; elements in a sequence are typed as biological
symbol instead of character or byte. They are special purpose types rather than
simply strings and hence offer additional functionality that naive string types
don't have. Though this strictness sacrifices some convenience, it also means
you can always rely on a DNA sequence type to store DNA and nothing but DNA,
without having to check, or deal with lowercase versus uppercase and so on.
Strict separation of sequence types also means we are free to choose the most
efficient representation. DNA and RNA sequences are encoded using either four
bits per base (which is the default), or two bits per base. This makes them
memory efficient and allows us to speed up many common operations and
transformations, like nucleotide composition, reverse complement, and *k*-mer
enumeration.

The `Bio.Seq` provides three different sequence types: `BioSequence`, `Kmer` and
`ReferenceSequence`. Each of these types is a subtype of an abstract type called
`Sequence` and supports various string-like operations such as random access and
iteration. Different sequence types have different features. In most situations,
`BioSequence` type will do and is used as the default representation. But
sometimes other types are much more preferable in terms of memory efficiency and
computation performance.  Here is the summary table of these three types:

| Type                       | Description                                | Element type          | Mutability  | Allocation       |
| :----                      | :-----------                               | :------------         | :---------- | :----------      |
| `BioSequence{A<:Alphabet}` | general-purpose biological sequences       | DNA, RNA, Amino acids | mutable     | heap             |
| `Kmer{T<:NucleicAcid,k}`   | specialized for short nucleotide sequences | DNA, RNA              | immutable   | stack / register |
| `ReferenceSequence`        | specialized for long reference genomes     | DNA                   | immutable   | heap             |

Details of these different representations are explained in the following
sections:

* `BioSequence`: [General-purpose sequences](@ref)
* `Kmer`: [Nucleic acid k-mers](@ref)
* `ReferenceSequence`: [Reference sequences](@ref)
