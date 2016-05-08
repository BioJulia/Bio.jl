# Tools: Misc tools and function wrappers

    {meta}
    CurrentModule = Bio.Blast
    
## Blast wrapper
The `Blast` module is a wrapper for the command line interface of [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) from NCBI. It requires that you have Blast+ [installed](https://www.ncbi.nlm.nih.gov/books/NBK279671/) and accessible in your PATH (eg. you should be able to execute `$ blastn -h` from the command line).

### The Basics

This module allows you to run protein and nucleotide BLAST (`blastp` and `blastn` respectively) within julia and to parse BLAST results into Bio.jl types.

```julia
using Bio.Seq,
      Bio.Blast

seq1 = dna"CGGACCAGACGGACACAGGGAGAAGCTAGTTTCTTTCATGTGATTGANATNATGACTCTACTCCTAAAAGGGAAAAANCAATATCCTTGTTTACAGAAGAGAAACAAACAAGCCCCACTCAGCTCAGTCACAGGAGAGAN"
seq2 = dna"CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATCTTACTCTTTCATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTTTACAGATGTGAAACTTTCAA"

blastn(seq1, seq2)
```

These functions return a `Vector{BlastResult}`. Each element is a hit which includes the sequence of the hit, an [`AlignedSequence`](http://biojulia.github.io/Bio.jl/latest/man/alignments/) using the original query as a reference and some additional information (expect vaue, bitscore) for the hit.

```julia
immutable BlastResult
    bitscore::Float64
    expect::Float64
    queryname::ASCIIString
    hitname::ASCIIString
    hit::BioSequence
    alignment::AlignedSequence
end
```

If you've already run a blast analysis or have downloaded blast results in XML format from NCBI you can also pass an XML string to `readblastXML()` in order to obtain an array of `BlastResult`s.

```julia
results = readall(open("blast_results.xml")) # need to use `readstring` instead of `readall` for v0.5
readblastXML(results)
```

When parsing protein blast results, you must include the argument `seqtype="prot"`, eg. `readblastXML("results, seqtype="prot")`

### Options for `blastn` and `blastp`

Both of the basic BLAST+ commands can accept a single `BioSequence`, a `Vector{BioSequence}` or a sting representing a file path to a fasta formatted file as arguments for both `query` and `subject`.

```julia
blastn([seq1, seq2], [seq2, seq3])

blastp(aaseq, "path/to/sequences.fasta")
```

If you have a local blast database (eg through the use of `$ makeblastdb`), you can use this database as the `subject`

```julia
blastn(seq1, "path/to/blast_db", db=true)
```

If you want to modify the search using additional options (eg. return only results with greater than 90% identity), you may pass a `Vector` of flags (see [here](http://www.ncbi.nlm.nih.gov/books/NBK279675/) for valid arguments - do not use flags that will alter file handling such as `-outfmt`)

```julia
blastn(seq1, seq2, ["-perc_identity", 90, "-evalue", "9.0"])
```
