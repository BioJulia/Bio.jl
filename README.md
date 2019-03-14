**NOTICE:**

**Since julia version 1.0 and the new Pkg system, the purpose of this package
has become moot. You can still download and use this package as a conveinience,
but consider instead managing the BioJulia packages you're using with `Pkg`
and its support for isolated projects. See the list below to find the
BioJulia packages that provide this package's modules now.**

**This package is inactive, meaning that updates to it will only be made as
and when time allows.**

---

# <img src="./sticker.svg" width="30%" align="right" /> Bio

[![Latest release](https://img.shields.io/github/release/BioJulia/Bio.jl.svg)](https://github.com/BioJulia/Bio.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/Bio.jl/blob/master/LICENSE)
[![Pkg Status](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)
[![Chat](https://img.shields.io/gitter/room/BioJulia/Bio.jl.svg)](https://gitter.im/BioJulia/Bio.jl)


## Description

Bio provides a framework for computational biology and bioinformatics.

This version of `Bio` package is actually better described as a meta-package.
It actually consolidates many other smaller packages in the BioJulia package
ecosystem and makes them easier to install and use together, with less worry
about version compatiblity and dependencies.

`Bio` has the current feature modules, go check out their documentation pages to
see how they work and how to use them:

* `Bio.Seq`: A biological sequences module (provided by [_BioSequences.jl_](https://github.com/BioJulia/BioSequences.jl))
    * Biological symbols (DNA, RNA and amino acid)
    * Biological sequences
    * Sequence search algorithms
    * Readers for FASTA, FASTQ and .2bit file formats
* `Bio.Align`: Sequence alignment (provided by [_BioAlignments.jl_](https://github.com/BioJulia/BioAlignments.jl))
    * Alignment data structures
    * Pairwise alignment algorithms
    * Reader for SAM and BAM file formats
* `Bio.Intervals`: Genomic intervals and annotations (provided by [_GenomicFeatures.jl_](https://github.com/BioJulia/GenomicFeatures.jl))
    * Genomic intervals with annotations
    * Readers for BED, bigWig, bigBed and GFF3 file formats
* `Bio.Structure`: Molecular structures (provided by [_BioStructures.jl_](https://github.com/BioJulia/BioStructures.jl))
    * Macromolecular structures (e.g. proteins)
    * Reader for PDB file format
* `Bio.Var`: Biological variation (provided by [_GeneticVariation.jl_](https://github.com/BioJulia/GeneticVariation.jl))
    * Mutation counting
    * Genetic and evolutionary distances
    * Readers for VCF and BCF file formats
* `Bio.Phylo`: Phylogenetics (provided by [_Phylogenies.jl_](https://github.com/BioJulia/Phylogenies.jl))
    * Phylogenetic trees
* `Bio.Services`: APIs to other services (provided by [_BioServices.jl_](https://github.com/BioJulia/BioServices.jl))
    * Entrez utilities (E-utilities)
* `Bio.Tools`: A module for running command line tools (provided by [_BioTools.jl_](https://github.com/BioJulia/BioTools.jl))
    * Run BLAST and parse its output


## Installation

Install Bio from the Julia REPL:

```julia
using Pkg
add("Bio")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Testing

Bio is a meta-package, and testing on CI infrastructure currently takes to long.
Each Bio release is tested on desktop before release, but if you run into problems,
please don't hesitate to contact a member of BioJulia.


## Contributing

We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing)
detailed contributor and maintainer guidelines, and code of conduct.


### Financial contributions

We also welcome financial contributions in full transparency on our
[open collective](https://opencollective.com/biojulia).
Anyone can file an expense. If the expense makes sense for the development
of the community, it will be "merged" in the ledger of our open collective by
the core contributors and the person who filed the expense will be reimbursed.


## Backers & Sponsors

Thank you to all our backers and sponsors!

Love our work and community? [Become a backer](https://opencollective.com/biojulia#backer).

[![backers](https://opencollective.com/biojulia/backers.svg?width=890)](https://opencollective.com/biojulia#backers)

Does your company use BioJulia? Help keep BioJulia feature rich and healthy by
[sponsoring the project](https://opencollective.com/biojulia#sponsor)
Your logo will show up here with a link to your website.

[![](https://opencollective.com/biojulia/sponsor/0/avatar.svg)](https://opencollective.com/biojulia/sponsor/0/website)
[![](https://opencollective.com/biojulia/sponsor/1/avatar.svg)](https://opencollective.com/biojulia/sponsor/1/website)
[![](https://opencollective.com/biojulia/sponsor/2/avatar.svg)](https://opencollective.com/biojulia/sponsor/2/website)
[![](https://opencollective.com/biojulia/sponsor/3/avatar.svg)](https://opencollective.com/biojulia/sponsor/3/website)
[![](https://opencollective.com/biojulia/sponsor/4/avatar.svg)](https://opencollective.com/biojulia/sponsor/4/website)
[![](https://opencollective.com/biojulia/sponsor/5/avatar.svg)](https://opencollective.com/biojulia/sponsor/5/website)
[![](https://opencollective.com/biojulia/sponsor/6/avatar.svg)](https://opencollective.com/biojulia/sponsor/6/website)
[![](https://opencollective.com/biojulia/sponsor/7/avatar.svg)](https://opencollective.com/biojulia/sponsor/7/website)
[![](https://opencollective.com/biojulia/sponsor/8/avatar.svg)](https://opencollective.com/biojulia/sponsor/8/website)
[![](https://opencollective.com/biojulia/sponsor/9/avatar.svg)](https://opencollective.com/biojulia/sponsor/9/website)


## Questions?

If you have a question about contributing or using BioJulia software, come
on over and chat to us on [Gitter](https://gitter.im/BioJulia/General), or you can try the
[Bio category of the Julia discourse site](https://discourse.julialang.org/c/domain/bio).