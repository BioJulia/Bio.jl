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

```@doc
Phylogeny
```

### Constructors

You can create a very simple unresolved phylogeny (a star phylogeny) by
providing the tips as a vector of strings or a vector of symbols.

```@example
tips = [:A, :B, :C]
tree = Phylogeny(tips)
```

```@example starstring
tips = ["A", "B", "C"]
tree = Phylogeny(tips)
```

### Roots

You can test whether a phylogeny is rooted, is re-rootable, and get the root
vertex of a phylogeny. You can also test if a vertex of a phylogeny is a root.

```@doc
isrooted
isrerootable
root
``

```@example
isrooted(tree)
```

```@example
isrerootable(tree)
```

```@example
root(tree)
```
