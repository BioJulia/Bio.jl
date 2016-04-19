# Bio.Phylo: Phylogenetic Trees

    {meta}
    CurrentModule = Bio.Phylo
    DocTestSetup = quote
        using Bio.Phylo
    end

The `Bio.Phylo` module provides several data types for handling phylogenetic
trees and models of ancestry/evolutionary history.

Phylogeneies are inferred models of the inferred evolutionary relationships of a
set of species, or strains, or populations. The inference is made usually based
on physical or genetic similarity.

Phylogenies can be rooted, or unrooted, and
phylogenies can also be vary varied as to the information they contain.
For example, a cladogram is a phylogeny formed using cladistic methods.
Cladograms describe the branching pattern of the history of a set of taxa.
Cladograms do not contain any branch lengths which represent the amount of time
passed, or amount of genetic change or evolutionary distance. To provide an
example a cladogram would tell us that we, branch closely with chimps, but will not
us how much evolutionary distance is between us, or how much time has passed since
the split.

A Phylogram is like a cladogram, but it _does_ include branch lengths that tell
you how much genetic change has occurred (typically as an evolutionary distance).

A Chronogram is like a Phylogram, but it's branch lengths indicate evolutionary
time.

Phylogenies may also have their nodes annotated with any number of pieces of extra
information and metadata.


## PhyNodes: The building blocks of phylogenies.
