# Bio.Phylo: Phylogenetic trees and networks

```@meta
CurrentModule = Bio.Phylo
DocTestSetup = quote
    using Bio.Phylo
end
```

The `Bio.Phylo` module is for data types and methods for handling phylogenetic
trees and networks.

## Phylogenies

```@docs
Phylogeny
```

### Constructors

You can create a very simple unresolved phylogeny (a star phylogeny) by
providing the tips as a vector of strings or a vector of symbols.

```@example phylo
using Bio.Phylo # hide
tips = [:A, :B, :C]
tree = Phylogeny(tips)
```

```@example
using Bio.Phylo # hide
tips = ["A", "B", "C"]
tree = Phylogeny(tips)
```

### Roots

You can test whether such a phylogeny is rooted, is re-rootable, and get the
root vertex of a phylogeny.
You can also test if a vertex of a phylogeny is a root.

```@docs
isrooted
isrerootable
root
```

```@example phylo
isrooted(tree)
```

```@example phylo
isrerootable(tree)
```

```@example phylo
root(tree)
```

## Divergence time estimation

`Phylo` has a submodule called `Dating` which contains methods for divergence
time estimation between sequences.

### Dating methods

Currently `Phylo.Dating` has two types which are used as function arguments to
dictate how to compute coalescence times. They all inherit from the abstract
data type `DatingMethod`.

```@docs
Dating.SimpleEstimate
Dating.SpeedDating
```
