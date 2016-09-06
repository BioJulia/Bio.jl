<p align="center"><img src="https://raw.githubusercontent.com/BioJulia/assets/master/branding/bio/BioJl_Design_1.png" width="50%" alt="Bio.jl" /></p>

As the flagship package of the BioJulia organisation, Bio.jl provides core
modules containing efficient data types and algorithms, that most
bioinformaticians and biologists would want to use for analyses or for building
their own applications. Bio.jl is built on top of the [Julia programming
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
* `Bio.Intervals`: Genomic intervals and annotations
    * Genomic intervals with annotations
    * Readers for BED and BigBed file formats
* `Bio.Structure`: Molecular structures
    * Macromolecular structures (e.g. proteins)
    * Reader for PDB file format
* `Bio.Var`: Biological variation
    * Mutation counting
    * Genetic and evolutionary distances
* `Bio.Phylo`: Phylogenetics
    * Phylogenetic trees


## Badges

Get Help: [![Join the chat at Gitter!](https://badges.gitter.im/BioJulia.png)](https://gitter.im/BioJulia/Bio.jl)
[![reference docs](https://img.shields.io/badge/docs-reference-blue.svg)](http://biojulia.github.io/Bio.jl/latest/)

Activity: [![Planned Work](https://badge.waffle.io/BioJulia/Bio.jl.svg?label=stage:%20planning&title=Planned)](http://waffle.io/BioJulia/Bio.jl)
[![Work In Progress](https://badge.waffle.io/BioJulia/Bio.jl.svg?label=stage:%20WIP&title=In%20Progress)](http://waffle.io/BioJulia/Bio.jl)

Code Quality: [![Build Status](https://travis-ci.org/BioJulia/Bio.jl.svg?branch=master)](https://travis-ci.org/BioJulia/Bio.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/nq4w264346py8esp/branch/master?svg=true)](https://ci.appveyor.com/project/Ward9250/bio-jl/branch/master)
[![Coverage Status](https://img.shields.io/coveralls/BioJulia/Bio.jl.svg)](https://coveralls.io/r/BioJulia/Bio.jl)
[![codecov.io](http://codecov.io/github/BioJulia/Bio.jl/coverage.svg?branch=master)](http://codecov.io/github/BioJulia/Bio.jl?branch=master)


## Documentation

Read the reference manual here: <http://biojulia.github.io/Bio.jl/latest>.


## Install

Bio.jl is a registered package in the official package management system and can
be installed with a command:
```julia
julia> Pkg.add("Bio")

```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Contributing

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features. If you have a question about
contributing, our [Gitter chat room](https://gitter.im/BioJulia/Bio.jl) would be
the best place to communicate with other users and developers.


## Roadmap

Our roadmap is on the wiki: <https://github.com/BioJulia/Bio.jl/wiki/roadmap>.
