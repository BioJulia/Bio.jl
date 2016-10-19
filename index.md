---
layout: pkgpage
title: Bio.jl
subtitle: Bioinformatics and Computational Biology infastructure for Julia
logo: img/Bio_Brand.svg
---

**Introduction**

As the flagship package of the BioJulia organisation, Bio.jl provides core
modules containing efficient data types and algorithms, that most
bioinformaticians and biologists would want to use for analyses or for building
their own applications.

Bio.jl is built on top of the Julia programming language, a high-level and
high-performance programming language for technical computing.
Bio.jl and Julia are open source and their source codes are immediately
available to the public.

---

**Overview**

_Modules and Functionality_

Bio.jl provides programmable components for quick prototyping of new analyses
and algorithms.
These components are carefully tuned to achieve the best performance without
sacrificing the usability of the dynamic programming language.
The following modules are currently part of the package and actively developed
as submodules:

<div class="row">

    <div class="col-md-4">
        <center>
            <img id="seqicon" src="img/module_icons/Biojl_Seq_Icon_Blue.svg" width="25%" />
            <br><strong>Seq</strong>
        </center>
        <ul>
            <li>Biological symbols (DNA, RNA, and amino acids)</li>
            <li>Biological sequences</li>
            <li>Sequence search algorithms</li>
            <li>Readers for FASTA, FASTQ and .2bit file formats</li>
        </ul>
    </div>

    <div class="col-md-4">
        <center>
            <img id="alignicon" src="img/module_icons/Biojl_Align_Icon_Green.svg" width="25%" />
            <br><strong>Align</strong>
        </center>
        <ul>
            <li>Biological sequence alignment</li>
            <li>Alignment data structures</li>
            <li>Pairwise alignment algorithms</li>
        </ul>
    </div>

    <div class="col-md-4">
        <center>
            <img id="intervalsicon" src="img/module_icons/Biojl_Intervals_Icon_Purple.svg" width="25%" />
            <br><strong>Intervals</strong>
        </center>
        <ul>
            <li>Intervals and annotations</li>
            <li>Genomic intervals with annotations</li>
            <li>Readers for BED and BigBed file formats</li>
        </ul>
    </div>

</div>

<hr>

<div class="row">

    <div class="col-md-4">
        <center>
            <img id="structicon" src="img/module_icons/Biojl_Struct_Icon_Red.svg" width="25%" />
            <br><strong>Structure</strong>
        </center>
        <ul>
            <li>Molecular structures</li>
            <li>Macromolecular structures (e.g. proteins)</li>
            <li>Reader for the PDB file format</li>
        </ul>
    </div>

    <div class="col-md-4">
        <center>
            <img id="varicon" src="img/module_icons/Biojl_Var_Icon_Blue.svg" width="25%" />
            <br><strong>Var</strong>
        </center>
        <ul>
            <li>Biological variation</li>
            <li>Mutation counting</li>
            <li>Genetic and evolutionary distances</li>
        </ul>
    </div>

    <div class="col-md-4">
        <center>
            <img id="phyloicon" src="img/module_icons/Biojl_Phylo_Icon_Green.svg" width="25%" />
            <br><strong>Phylo</strong>
        </center>
        <ul>
            <li>Phylogenetic trees</li>
        </ul>
    </div>

</div>

---

**Quick Start**

In a Julia console, enter:

{% highlight julia %}
Pkg.add("Bio.jl")
{% endhighlight %}

Julia's package manager will try to find the latest released version and
install it.

Alternatively, to install the bleeding edge of a branch use:

{% highlight julia %}
Pkg.clone("https://github.com/BioJulia/Bio.jl.git")
{% endhighlight %}

Then use the Bio submodules in a script like so:

{% highlight julia %}
using Bio.Seq
{% endhighlight %}

<script src="js/biojl.js"></script>
