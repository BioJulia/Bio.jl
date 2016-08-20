# Bio.Var: Biological Variation.

```@meta
CurrentModule = Bio.Var
DocTestSetup = quote
    using Bio.Var
end
```

## Counting mutations

You can count the numbers of different types of mutations between two aligned
(i.e. of equal length) nucleotide sequences.

### Different types of mutation

The types of mutations that can currently be counted are displayed in the
following table.

| Type                   | Meaning                                 |
| :--------------------- | :-------------------------------------- |
| `DifferentMutation`    | Nucleotides that are different          |
| `TransitionMutation`   | Transition mutations                    |
| `TransversionMutation` | Transversion mutations                  |

`DifferentMutation` simply compares two nucleotides to see if they can be deemed
the same or different. Ambiguous cases are detected and excluded from the
computation.

`TransitionMutation` compares two nucleotides to see if they can be determined
to be a transition mutation.
As with `DifferentMutation`, cases where a nucleotide is ambiguous are excluded
from the computation.

`TransversionMutation` compares two nucleotides to see if they can be determined
to be a transversion mutation.
As with `DifferentMutation`, cases where a nucleotide is ambiguous are excluded
from the computation.

### `count_mutations` method

Mutations are counted using the `count_mutations` method.
The method outputs a tuple. The first value is the number of mutations counted.
The second value is the number of sites examined. Sites which have gaps and
uncertain nucleotides are not examined and so this second value will be less
than the length of the two biological sequences.

```julia
count_mutations(dna"ATCGATCG", dna"ACCGATCG", DifferentMutation)

count_mutations(dna"ATCGATCG", dna"ACCGATCG", TransitionMutation)

count_mutations(dna"ATCGATCG", dna"ACCGATCG", TransversionMutation)

count_mutations(dna"ATCGATCG", dna"ACCGATCG", TransitionMutation, TransversionMutation)

count_mutations(rna"AUCGAUCG", rna"ACCGAUCG", DifferentMutation)

count_mutations(rna"AUCGAUCG", rna"ACCGAUCG", TransitionMutation)

count_mutations(rna"AUCGAUCG", rna"ACCGAUCG", TransversionMutation)

count_mutations(rna"AUCGAUCG", rna"ACCGAUCG", TransitionMutation, TransversionMutation)
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

Using the distance method you can compute different kinds of evolutionary
sequence.

```@docs
distance(Count{T}, BioSequence, BioSequence)
distance(Proportion{T}, BioSequence, BioSequence)
distance(JukesCantor69, BioSequence, BioSequence)
distance(Kimura80, BioSequence, BioSequence)
```
