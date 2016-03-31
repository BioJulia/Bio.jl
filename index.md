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

---

**Overview**

_Modules and Functionality_

The following modules and functionality are currently part of Bio.jl and/or are
being actively developed, as submodules of the Bio module.

- Biological sequences: _Seq_

- Biological sequence alignment: _Align_

- Intervals and annotations: _Intervals_

- Molecular structures: _Structure_

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
