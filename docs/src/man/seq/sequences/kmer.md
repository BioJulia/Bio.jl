```@meta
CurrentModule = Bio.Seq
DocTestSetup = quote
    using Bio.Seq
end
```

## Nucleic acid k-mers

A common strategy to simplify the analysis of sequence data is to operate or
short k-mers, for size fixed size `k`. These can be packed into machine integers
allowing extremely efficient code. The `Bio.Seq` module has built in support for
representing short sequences in 64-bit integers. Besides being fixed length,
`Kmer` types, unlike other sequence types cannot contain ambiguous symbols like
'N'.

The `Kmer{T,k}` type parameterized on symbol type (`T`, either `DNA`,
or `RNA`) and size `k`. For ease of writing code, two type aliases for
each nucleotide type are defined and named as `DNAKmer{k}` and `RNAKmer{k}`:
```jldoctest
julia> DNAKmer("ACGT")  # create a DNA 4-mer from a string
DNA 4-mer:
ACGT

julia> RNAKmer("ACGU")  # create an RNA 4-mer from a string
RNA 4-mer:
ACGU

julia> kmer"ACGT" # DNA k-mers may also be written as literals
DNA 4-mer:
ACGT

julia> typeof(DNAKmer("ACGT"))
Bio.Seq.Kmer{BioSymbols.DNA,4}
```


```@docs
each
canonical
neighbors
```
