# Bio.Var: Biological Variation.

```@meta
CurrentModule = Bio.Var
DocTestSetup = quote
    using Bio.Var
end
```

## Counting mutations

You can count the numbers of different types of mutations in a pairwise manner
for a set of nucleotide sequences.

### Different types of mutation

The types of mutations that can currently be counted are `DifferentMutation`,s
`TransitionMutation`s, and `TransversionMutation`s.

```@docs
DifferentMutation
TransitionMutation
TransversionMutation
```

### `count_mutations` method

Mutations are counted using the `count_mutations` method.
The method outputs a tuple. The first value is the number of mutations counted.
The second value is the number of sites examined. Sites which have gaps and
uncertain nucleotides are not examined and so this second value will be less
than the length of the two biological sequences.

```julia
count_mutations(DifferentMutation, [dna"ATCGATCG", dna"ACCGATCG"])

count_mutations(TransitionMutation, [dna"ATCGATCG", dna"ACCGATCG"])

count_mutations(TransversionMutation, [dna"ATCGATCG", dna"ACCGATCG"])

count_mutations(TransitionMutation, TransversionMutation, [dna"ATCGATCG", dna"ACCGATCG"])

count_mutations(DifferentMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])

count_mutations(TransitionMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])

count_mutations(TransversionMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])

count_mutations(TransitionMutation, TransversionMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])
```

## Computing evolutionary and genetic distances

Just as you can count the number of mutations between two nucleotide sequences,
you can compute the evolutionary distance between two nucleotide sequences.

### Different evolutionary distance measures

The types of distances that can currently be computed are described below.

```@docs
Count{T}
Proportion{T}
JukesCantor69
Kimura80
```

### The distance method

```@docs
distance
```
