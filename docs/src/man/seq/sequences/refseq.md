## Reference sequences

`DNASequence` (alias of `BioSequence{DNAAlphabet{4}}`) is a flexible data
structure but always consumes 4 bits per base, which will waste a large part of
the memory space when storing reference genome sequences.  In such a case,
`ReferenceSequence` is helpful because it compresses positions of 'N' symbols so
that long DNA sequences are stored with almost 2 bits per base. An important
limitation is that the `ReferenceSequence` type is immutable due to the
compression. Other sequence-like operations are supported:
```jlcon
julia> seq = ReferenceSequence(dna"NNCGTATTTTCN")
12nt Reference Sequence:
NNCGTATTTTCN

julia> seq[1]
DNA_N

julia> seq[5]
DNA_T

julia> seq[2:6]
5nt Reference Sequence:
NCGTA

julia> ReferenceSequence(dna"ATGM")  # DNA_M is not accepted
ERROR: ArgumentError: invalid symbol M ∉ {A,C,G,T,N} at 4
 in convert at /Users/kenta/.julia/v0.4/Bio/src/seq/refseq.jl:58
 in call at essentials.jl:56

```
