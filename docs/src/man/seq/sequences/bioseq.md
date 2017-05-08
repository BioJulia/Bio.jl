```@meta
CurrentModule = Bio.Seq
DocTestSetup = quote
    using Bio.Seq
end
```
## General-purpose sequences

`BioSequence{A}` is a generic sequence type parameterized by an alphabet type
`A` that defines the domain (or set) of biological symbols, and each alphabet
has an associated symbol type. For example, `AminoAcidAlphabet` is associated
with `AminoAcid` and hence an object of the `BioSequence{AminoAcidAlphabet}`
type represents a sequence of amino acids.  Symbols from multiple alphabets
can't be intermixed in one sequence type.

The following table summarizes common sequence types that are defined in the
`Bio.Seq` module:

| Type                               | Symbol type | Type alias          |
| :--------------------------------- | :---------- | :------------------ |
| `BioSequence{DNAAlphabet{4}}`      | `DNA`       | `DNASequence`       |
| `BioSequence{RNAAlphabet{4}}`      | `RNA`       | `RNASequence`       |
| `BioSequence{AminoAcidAlphabet}`   | `AminoAcid` | `AminoAcidSequence` |
| `BioSequence{CharAlphabet}`        | `Char`      | `CharSequence`      |

Parameterized definition of the `BioSequence{A}` type is for the purpose of
unifying the data structure and operations of any symbol type. In most cases,
users don't have to care about it and can use *type aliases* listed above.
However, the alphabet type fixes the internal memory encoding and plays an
important role when optimizing performance of a program
(see [Using a more compact sequence representation](@ref) section for low-memory
encodings).  It also enables a user to define their own alphabet only by
defining few numbers of methods.
This is described in [Defining a new alphabet](@ref) section.


### Constructing sequences

#### Using string literals

Most immediately, sequence literals can be constructed using the string macros
`dna`, `rna`, `aa`, and `char`:

```jldoctest
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

However it should be noted that by default these sequence literals
allocate the `BioSequence` object before the code containing the sequence
literal is run.
This means there may be occasions where your program does not behave as you
first expect.
For example consider the following code:

```jldoctest
julia> function foo()
           s = dna"CTT"
           push!(s, DNA_A)
       end
foo (generic function with 1 method)

```

```@meta
DocTestSetup = quote
    using Bio.Seq
    function foo()
        s = dna"CTT"d
        push!(s, DNA_A)
    end
end
```

You might expect that every time you call `foo`, that a DNA sequence `CTTA` would
be returned. You might expect that this is because every time `foo` is called,
a new DNA sequence variable `CTT` is created, and and `A` nucleotide is pushed
to it, and the result, `CTTA` is returned.
In other words you might expect the following output:

```jldoctest
julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

```

However, this is not what happens, instead the following happens:

```@meta
DocTestSetup = quote
    using Bio.Seq
    function foo()
        s = dna"CTT"s
        push!(s, DNA_A)
    end
end
```

```jldoctest
julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
5nt DNA Sequence:
CTTAA

julia> foo()
6nt DNA Sequence:
CTTAAA

```

The reason for this is because the sequence literal is allocated only once
before the first time the function `foo` is called and run. Therefore, `s` in
`foo` is always a reference to that one sequence that was allocated.
So one sequence is created before `foo` is called, and then it is pushed to
every time `foo` is called. Thus, that one allocated sequence grows with every
call of `foo`.

If you wanted `foo` to create a new sequence each time it is called,
then you can add a flag to the end of the sequence literal to dictate behaviour:
A flag of 's' means 'static': the sequence will be allocated before code is run,
as is the default behaviour described above.
However providing 'd' flag changes the behaviour: 'd' means 'dynamic':
the sequence will be allocated at whilst the code is running, and not before.
So to change `foo` so as it creates a new sequence
each time it is called, simply add the 'd' flag to the sequence literal:
```@meta
DocTestSetup = quote
    using Bio.Seq
end
```

```jldoctest
julia> function foo()
           s = dna"CTT"d     # 'd' flag appended to the string literal.
           push!(s, DNA_A)
       end
foo (generic function with 1 method)

```

Now every time `foo` is called, a new sequence `CTT` is created, and an `A`
nucleotide is pushed to it:

```@meta
DocTestSetup = quote
    using Bio.Seq
    function foo()
        s = dna"CTT"d
        push!(s, DNA_A)
    end
end
```

```jldoctest
julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

julia> foo()
4nt DNA Sequence:
CTTA

```

So the take home message of sequence literals is this:

Be careful when you are using sequence literals inside of functions, and inside
the bodies of things like for loops. And if you use them and are unsure, use the
 's' and 'd' flags to ensure the behaviour you get is the behaviour you intend.


#### Other constructors and conversion

Sequences can also be constructed from strings or arrays of nucleotide or amino
acid symbols using constructors or the `convert` function:

```jldoctest
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
```jldoctest
julia> convert(String, dna"TTANGTA")
"TTANGTA"

julia> convert(Vector{DNA}, dna"TTANGTA")
7-element Array{BioSymbols.DNA,1}:
 DNA_T
 DNA_T
 DNA_A
 DNA_N
 DNA_G
 DNA_T
 DNA_A

```

Sequences can also be concatenated into longer sequences:
```jldoctest
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
```jldoctest
julia> dna = dna"TTANGTAGACCG"
12nt DNA Sequence:
TTANGTAGACCG

julia> rna = convert(RNASequence, dna)
12nt RNA Sequence:
UUANGUAGACCG

julia> dna.data === rna.data  # underlying data are same
true

```

A random sequence can be obtained by the `randdnaseq`, `randrnaseq` and
`randaaseq` functions, which generate `DNASequence`, `RNASequence` and
`AminoAcidSequence`, respectively. Generated sequences are composed of the
standard symbols without ambiguity and gap. For example, `randdnaseq(6)` may
generate `dna"TCATAG"` but never generates `dna"TNANAG"` or `dna"T-ATAG"`.

A translatable `RNASequence` can also be converted to an `AminoAcidSequence`
using the [`translate`](@ref) function.


### Indexing, modifying and transformations

#### Getindex

Sequences for the most part behave like other vector or string types. They can
be indexed using integers or ranges:

```jldoctest
julia> seq = dna"ACGTTTANAGTNNAGTACC"
19nt DNA Sequence:
ACGTTTANAGTNNAGTACC

julia> seq[5]
DNA_T

julia> seq[6:end]
14nt DNA Sequence:
TANAGTNNAGTACC

```

Note that, indexing a biological sequence by range creates a subsequence of the
original sequence. Unlike `Arrays` in the standard library, creating a
subsequence is copy-free: a subsequence simply points to the original sequence
data with its range. You may think that this is unsafe because modifying
subsequences propagates to the original sequence, but this doesn't happen
actually:

```jldoctest
julia> seq = dna"AAAA"    # create a sequence
4nt DNA Sequence:
AAAA

julia> subseq = seq[1:2]  # create a subsequence from `seq`
2nt DNA Sequence:
AA

julia> subseq[2] = DNA_T  # modify the second element of it
DNA_T

julia> subseq             # the subsequence is modified
2nt DNA Sequence:
AT

julia> seq                # but the original sequence is not
4nt DNA Sequence:
AAAA

```

This is because modifying a sequence checks whether its underlying data are
shared with other sequences under the hood. If and only if the data are shared,
the subsequence creates a copy of itself. Any modifying operation does this
check. This is called *copy-on-write* strategy and users don't need to care
about it because it is transparent: If the user modifies a sequence with or
subsequence, the job of managing and protecting the underlying data of sequences
is handled for them.


#### Setindex and modifying DNA sequences

The biological symbol at a given locus in a biological sequence can be set using
setindex:

```jldoctest
julia> seq = dna"ACGTTTANAGTNNAGTACC"
19nt DNA Sequence:
ACGTTTANAGTNNAGTACC

julia> seq[5] = DNA_A
DNA_A

```

In addition, many other modifying operations are possible for biological
sequences such as `push!`, `pop!`, and `insert!`, which should be familiar to
people used to editing arrays.

```@docs
push!
pop!
shift!
unshift!
insert!
deleteat!(::Bio.Seq.BioSequence, ::Integer)
append!
copy!
```

Here are some examples:

```jldoctest
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

julia> deleteat!(seq, 2)
5nt DNA Sequence:
AGTAT

julia> deleteat!(seq, 2:3)
3nt DNA Sequence:
AAT

```

#### Additional transformations

In addition to these basic modifying functions, other sequence transformations
which are common in bioinformatics are also provided.

```@docs
reverse!
complement!
reverse_complement!
```

```jldoctest
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

##### Translation

Translation is a slightly more complex transformation for RNA Sequences and so
we describe it here in more detail.

The [`translate`](@ref) funtion translates a sequence of codons in a RNA sequence
to a amino acid sequence besed on a genetic code mapping. The `Bio.Seq` module
contains all NCBI defined genetic codes and they are registered in
[`ncbi_trans_table`](@ref).

```@docs
translate
ncbi_trans_table
```

```jldoctest
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


### Iteration

Sequences also work as iterators over symbols:

```jldoctest
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


### Using a more compact sequence representation

As we saw above, DNA and RNA sequences can store any ambiguous nucleotides like
'N'.  If you are sure that nucleotide sequences store unambiguous nucleotides
only, you can save the memory space of sequences. `DNAAlphabet{2}` is an
alphabet that uses two bits per base and limits to only unambiguous nucleotide
symbols (ACGT in DNA and ACGU in RNA). To create a sequence of this
alphabet, you need to explicitly pass `DNAAlphabet{2}` to `BioSequence` as its
type parameter:

```jldoctest
julia> seq = BioSequence{DNAAlphabet{2}}("ACGT")
4nt DNA Sequence:
ACGT

```

Recall that `DNASequence` is a type alias of `BioSequence{DNAAlphabet{4}}`,
which uses four bits per base. That is, `BioSequence{DNAAlphabet{2}}` saves half
of memory footprint compared to `BioSequence{DNAAlphabet{4}}`. If you need to
handle reference genomes that are composed of five nucleotides, ACGTN,
consider to use the `ReferenceSequence` type described in the [Reference
sequences](@ref) section.


### Defining a new alphabet

The alphabet type parameter `A` of `BioSequence{A}` enables a user to extend
functionality of `BioSequence` with minimum effort. As an example, definition of
a new alphabet type representing a sequence of boolean values is shown below:

```jldoctest
julia> immutable BoolAlphabet <: Alphabet end

julia> Seq.bitsof(::Type{BoolAlphabet}) = 1

julia> Seq.eltype(::Type{BoolAlphabet}) = Bool

julia> Seq.alphabet(::Type{BoolAlphabet}) = false:true

julia> function Seq.encode(::Type{BoolAlphabet}, x::Bool)
           return UInt64(ifelse(x, 0x01, 0x00))
       end

julia> function Seq.decode(::Type{BoolAlphabet}, x::UInt64)
           if x > 0x01
               throw(Seq.DecodeError(BoolAlphabet, x))
           end
           return ifelse(x == 0x00, false, true)
       end

```
