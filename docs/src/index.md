# Bio.jl Documentation

This documentation site contains [The Bio.jl User Manual](@ref)

As the flagship package of the BioJulia organisation, Bio.jl provides core
modules containing efficient data types and algorithms, that most
bioinformaticians and biologists would want to use for analyses or for building
their own applications.
Bio.jl is built on top of the Julia programming language, a high-level and
high-performance programming language for technical computing. Bio.jl and Julia
are open source and their source codes are immediately available to the public.

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
    * Readers for BED, BigBed and GFF3 file formats
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
