# Bio.Phylo: Phylogenetic trees and models of ancestry

```@meta
CurrentModule = Bio.Phylo
DocTestSetup = quote
    using Bio.Phylo
end
```

The `Bio.Phylo` module provides data types and methods for working with
phylogenetic trees and other similar structures that model evolutionary history
and relatedness between taxa.

## Types

### Phylogenies

The phylogeny datatype represents phylogenetic trees in BioJulia.

The internal structure and representation of the Phylogeny type uses a DiGraph.
DiGraph is a type representing a directed graph and is implemented in the
LightGraphs.jl package.
When a phylogenetic tree is created, the DiGraph is initiated with enough vertices
to be able to represent a fully bifurcating phylogeny, no more and no less.
The vertices of the DiGraph are reserved for different specific purposes.

Vertex ranges:
* 1:tree.ntaxa => leaf vertices
* tree.ntaxa + 1 => root vertex
* tree.ntaxa + 1:<#nodes> => clade vertices

As it is common to annotate phylogenies with data, the type is parametric with
two parameters C, and B. Parameter C represents a datatype that annotates clades
of the phylogeny, and parameter B represents a datatype that annotates branches
of the phylogeny.

Phylogenies can be created with the provided constructors for example,
constructing a phylogeny from just an array of Taxa names creates a basic and
fully unresolved star phylogeny, which can annotate clades and branches with
`Float64` values, which enables you to annotate branch lengths to branches and
confidence values to clades.

```jlcon
julia> tree = Phylogeny([:Human, :Chimp, :Dog, :Fish])
Bio.Phylo.Phylogeny{Float64,Float64}({7, 4} directed graph,Bio.Indexers.Indexer{Int64}(Dict(:Fish=>4,:Dog=>3,:Chimp=>2,:Root=>5,:Human=>1),[:Human,:Chimp,:Dog,:Fish,:Root]),Float64[],Dict{Pair{Int64,Int64},Float64}(),4,false,true)

julia> tree = Phylogeny(["Human", "Chimp", "Dog", "Fish"])
Bio.Phylo.Phylogeny{Float64,Float64}({7, 4} directed graph,Bio.Indexers.Indexer{Int64}(Dict(:Fish=>4,:Dog=>3,:Chimp=>2,:Root=>5,:Human=>1),[:Human,:Chimp,:Dog,:Fish,:Root]),Float64[],Dict{Pair{Int64,Int64},Float64}(),4,false,true)
```
