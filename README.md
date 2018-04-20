<p align="center"><img src="https://raw.githubusercontent.com/BioJulia/assets/master/branding/bio/optimised/BioJl_Design_1.png" width="50%" alt="Bio.jl" /></p>

_The Bioinformatics and Computation Biology infrastructure for the Julia language._

# The meta-package Bio.jl is currently undergoing restructuring and is not recommended for julia v0.6 users at this time, ask us on [![][gitter-img]][gitter-url] for more info, or check out this [PR](https://github.com/BioJulia/Bio.jl/issues/425). Individual BioJulia packages are stable on julia v0.6, see http://biojulia.net.


| **Chat** | **Documentation** | **Build Status** |
|:--------:|:-----------------:|:----------------:|
| [![][gitter-img]][gitter-url] | [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |


## Description

As the flagship package of the BioJulia organisation, Bio.jl collects core
modules containing efficient data types and algorithms that most
bioinformaticians and biologists would want to use for analyses or for building
their own applications.
These modules are generally from individual packages in the BioJulia ecosystem,
and Bio.jl is a convenient way to access them all in one place.
Bio.jl is built on top of the [Julia programming
language](http://julialang.org/), a high-level and high-performance programming
language for technical computing. Bio.jl and Julia are open source and their
source codes are immediately available to the public.

Bio.jl provides programmable components for quick prototyping of new analyses
and algorithms. These components are carefully tuned to achieve the best
performance without sacrificing the usability of the dynamic programming
language. The following modules are currently part of the package and actively
developed as submodules:
* `Bio.Seq`: Biological sequences
    * Biological symbols (DNA, RNA and amino acid)
    * Biological sequences
    * Sequence search algorithms
    * Readers for FASTA, FASTQ and .2bit file formats
* `Bio.Align`: Sequence alignment
    * Alignment data structures
    * Pairwise alignment algorithms
    * Reader for SAM and BAM file formats
* `Bio.Intervals`: Genomic intervals and annotations
    * Genomic intervals with annotations
    * Readers for BED, bigWig, bigBed and GFF3 file formats
* `Bio.Structure`: Molecular structures
    * Macromolecular structures (e.g. proteins)
    * Reader for PDB file format
* `Bio.Var`: Biological variation
    * Mutation counting
    * Genetic and evolutionary distances
    * Readers for VCF and BCF file formats
* `Bio.Phylo`: Phylogenetics
    * Phylogenetic trees
* `Bio.Services`: APIs to other services
    * Entrez utilities (E-utilities)


## Installation

Bio.jl is a registered package in the official package management system and can
be installed with a command:
```julia
julia> Pkg.add("Bio")

```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Documentation

- [**STABLE**][docs-stable-url] &mdash; **Documentation of the most recent stable release of Bio.jl.**
- [**LATEST**][docs-latest-url] &mdash; *Documentation of the in-development state of Bio.jl.*


## Contributing and Questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.

If you have a question about
contributing or using this package, our [Gitter chat room][gitter-url] would be
the best starting place to communicate with other users and developers.
You are encouraged to use the [Bio category of the Julia discourse
site](https://discourse.julialang.org/c/domain/bio) for technical questions.


## Roadmap

Our roadmap is on the wiki: <https://github.com/BioJulia/Bio.jl/wiki/roadmap>.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://biojulia.github.io/Bio.jl/latest
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://biojulia.github.io/Bio.jl/stable

[pkg-0.4-img]: http://pkg.julialang.org/badges/Bio_0.4.svg
[pkg-0.4-url]: http://pkg.julialang.org/?pkg=Bio
[pkg-0.5-img]: http://pkg.julialang.org/badges/Bio_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=Bio
[pkg-0.6-img]: http://pkg.julialang.org/badges/Bio_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=Bio

[codecov-img]: https://codecov.io/gh/BioJulia/Bio.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/BioJulia/Bio.jl
[travis-img]: https://travis-ci.org/BioJulia/Bio.jl.svg?branch=master
[travis-url]: https://travis-ci.org/BioJulia/Bio.jl
[appveyor-img]: https://ci.appveyor.com/api/projects/status/nq4w264346py8esp/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/Ward9250/bio-jl/branch/master

[gitter-img]: https://badges.gitter.im/BioJulia.png
[gitter-url]: https://gitter.im/BioJulia/Bio.jl
