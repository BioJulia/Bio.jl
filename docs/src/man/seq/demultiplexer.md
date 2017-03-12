
```@meta
CurrentModule = Bio.Seq
DocTestSetup = quote
    using Bio.Seq
end
```

# Sequence demultiplexing

Multiplex sequencing is a technology to sequence multiple samples at the same
time on a high-throughput DNA sequencer. Samples are distinguished by the short
prefix of a DNA sequence called DNA barcode. The `Bio.Seq` offers the
`Demultiplexer` type and the `demultiplex` function to identify the DNA barcode
of a longer DNA sequence allowing small errors.

In the following example, four kinds of DNA sequences of length 4 are used as
DNA barcodes. `Demultiplexer` takes these barcodes as its first argument with
a few options:
```jldoctest
julia> barcodes = DNASequence["ATGG", "CAGA", "GGAA", "TACG"];

julia> dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming)
Bio.Seq.Demultiplexer{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}}}:
  distance: hamming
  number of barcodes: 4
  number of correctable errors: 1
```

```@meta
DocTestSetup = quote
    using Bio.Seq
    barcodes = DNASequence["ATGG", "CAGA", "GGAA", "TACG"];
    dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming);
end
```

`n_max_errors` specifies the number of maximum correctable errors in a barcode.
The type of correctable errors depends on the `distance` parameter. When
`distance = :hamming` as shown above only substitutions are correctable. When
`distance = :levenshtein` substitutions, deletions, and insertions are
correctable. The user is responsible for keeping enough distances among
barcodes; `Demultiplexer` will throw an exception if two barcodes are within
`n_max_errors * 2`.

The `demultiplex` function takes a demultiplexer object and a DNA sequence, and
returns a tuple of a barcode index and a distance between the original barcode
sequence and the prefix sequence:
```jldoctest
julia> demultiplex(dplxr, dna"ATGGCGNT")  # 1st barcode with no errors
(1,0)

julia> demultiplex(dplxr, dna"CAGGCGNT")  # 2nd barcode with one error
(2,1)

julia> demultiplex(dplxr, dna"GGAACGNT")  # 3rd barcode with no errors
(3,0)

julia> demultiplex(dplxr, dna"TGACCGNT")  # no matching barcode
(0,-1)

```

The optional third argument controls the search strategy. `demultiplex` uses an
index to search the closest barcode within `n_max_errors` in the barcode set and
returns it if any by default. If the third argument is `true` it falls back to a
linear search after the index search and returns one of the closest barcodes at
random. The next example shows the difference of these two strategies:
```jldoctest
julia> demultiplex(dplxr, dna"TGACCGNT", false)  # linear search off (default)
(0,-1)

julia> demultiplex(dplxr, dna"TGACCGNT", true)   # linear search on
(3,2)

```
