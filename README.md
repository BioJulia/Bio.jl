# Bio.jl

| **Release**                                                     | **Documentation**                                                               | **Maintainers**                             |
|:---------------------------------------------------------------:|:-------------------------------------------------------------------------------:|:-------------------------------------------:|
| [![](https://img.shields.io/github/release/BioJulia/Bio.jl.svg)](https://github.com/BioJulia/Bio.jl/releases/latest) [![](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/Bio.jl/blob/master/LICENSE) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/Bio.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://biojulia.github.io/Bio.jl/latest) | ![](https://img.shields.io/badge/BioJulia%20Maintainer-Ward9250-orange.svg) |


## Description

| **Chat** | **Documentation** | **Build Status** | **Support** |
|:--------:|:-----------------:|:----------------:|:----------------:|
| [![][gitter-img]][gitter-url] | [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] | [![Backers on Open Collective](oc-backers-img)](#backers)[![Sponsors on Open Collective](oc-sponsors-img)](#sponsors) |
Bio.jl is the flagship meta-package of the BioJulia group.

Bio.jl re-exports and makes easy to install a set of BioJulia's most important
core bioinformatics packages.

The following modules are currently part of the meta-package:

* `Bio.Seq`: Biological sequences (provided by _BioSequences.jl_)
    * Biological symbols (DNA, RNA and amino acid)
    * Biological sequences
    * Sequence search algorithms
    * Readers for FASTA, FASTQ and .2bit file formats
* `Bio.Align`: Sequence alignment (provided by _BioAlignments.jl_)
    * Alignment data structures
    * Pairwise alignment algorithms
    * Reader for SAM and BAM file formats
* `Bio.Intervals`: Genomic intervals and annotations (provided by _GenomicFeatures.jl_)
    * Genomic intervals with annotations
    * Readers for BED, bigWig, bigBed and GFF3 file formats
* `Bio.Structure`: Molecular structures (provided by _BioStructures.jl_)
    * Macromolecular structures (e.g. proteins)
    * Reader for PDB file format
* `Bio.Var`: Biological variation (provided by _GeneticVariation.jl_)
    * Mutation counting
    * Genetic and evolutionary distances
    * Readers for VCF and BCF file formats
* `Bio.Phylo`: Phylogenetics (provided by _Phylogenies.jl_)
    * Phylogenetic trees
* `Bio.Services`: APIs to other services (provided by _BioServices.jl_)
    * Entrez utilities (E-utilities)


## Installation

Install Bio.jl from the Julia REPL:

```julia
using Pkg
add("Bio")
# Pkg.add("Bio") on julia v0.6-
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Testing

Bio.jl is tested against Julia `0.6` and current `0.7-dev` on Linux, OS X, and Windows.

| **PackageEvaluator**                                            | **Latest Build Status**                                                                                |
|:---------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------:|
| [![](https://pkg.julialang.org/badges/Bio_0.6.svg)](https://pkg.julialang.org/detail/Bio) [![](https://pkg.julialang.org/badges/Bio_0.7.svg)](https://pkg.julialang.org/detail/Bio) | [![](https://img.shields.io/travis/BioJulia/Bio.jl/master.svg?label=Linux+/+macOS)](https://travis-ci.org/BioJulia/Bio.jl) [![](https://ci.appveyor.com/api/projects/status/nq4w264346py8esp?svg=true)](https://ci.appveyor.com/project/Ward9250/bio-jl/branch/master) [![](https://codecov.io/gh/BioJulia/Bio.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/BioJulia/Bio.jl) |


## Contributing and Questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.
Please go to the [contributing section of the documentation](https://biojulia.net/Contributing/latest)
for more information.

If you have a question about
contributing or using this package, our [Gitter chat room][gitter-url] would be
the best starting place to communicate with other users and developers.
You are encouraged to use the [Bio category of the Julia discourse
site](https://discourse.julialang.org/c/domain/bio) for technical questions.


## Roadmap

Our roadmap is on the wiki: <https://github.com/BioJulia/Bio.jl/wiki/roadmap>.


## Backers

Thank you to all our backers! 🙏 [[Become a backer](https://opencollective.com/biojulia#backer)]

<a href="https://opencollective.com/biojulia#backers" target="_blank"><img src="https://opencollective.com/biojulia/backers.svg?width=890"></a>


## Sponsors

Support this project by becoming a sponsor. Your logo will show up here with a link to your website. [[Become a sponsor](https://opencollective.com/biojulia#sponsor)]

<a href="https://opencollective.com/biojulia/sponsor/0/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/0/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/1/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/1/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/2/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/2/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/3/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/3/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/4/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/4/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/5/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/5/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/6/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/6/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/7/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/7/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/8/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/8/avatar.svg"></a>
<a href="https://opencollective.com/biojulia/sponsor/9/website" target="_blank"><img src="https://opencollective.com/biojulia/sponsor/9/avatar.svg"></a>

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

[oc-backers-img]: https://opencollective.com/biojulia/backers/badge.svg
[oc-sponsors-img]:https://opencollective.com/biojulia/sponsors/badge.svg