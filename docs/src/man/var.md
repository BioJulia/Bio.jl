# Bio.Var: Biological variance module.

```@meta
CurrentModule = Bio.Var
DocTestSetup = quote
    using Bio.Var
end
```

### Counting mutations

You can count the numbers of different types of mutations between two aligned
DNA sequences.

#### Different types of mutation

The types of mutations that can be counted are currently.

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

#### `count_mutations` method

Mutations are counted using the `count_mutations` method.
