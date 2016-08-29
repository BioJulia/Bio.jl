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
count_mutations([dna"ATCGATCG", dna"ACCGATCG"], DifferentMutation)

count_mutations([dna"ATCGATCG", dna"ACCGATCG"], TransitionMutation)

count_mutations([dna"ATCGATCG", dna"ACCGATCG"], TransversionMutation)

count_mutations([dna"ATCGATCG", dna"ACCGATCG"], TransitionMutation, TransversionMutation)

count_mutations([rna"AUCGAUCG", rna"ACCGAUCG"], DifferentMutation)

count_mutations([rna"AUCGAUCG", rna"ACCGAUCG"], TransitionMutation)

count_mutations([rna"AUCGAUCG", rna"ACCGAUCG"], TransversionMutation)

count_mutations([rna"AUCGAUCG", rna"ACCGAUCG"], TransitionMutation, TransversionMutation)
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
