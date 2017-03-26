# Bio.Var: Biological Variation.

```@meta
CurrentModule = Bio.Var
DocTestSetup = quote
    using Bio.Var
end
```

## Identifying and counting sequence sites

You can identify and count the number of various types of site in a nucleotide
sequence, or pair/set of aligned nucleotide sequences.

### Different types of site

#### The abstract Site types

```@docs
Site
Mutation
```

#### The concrete types of site case

The types of mutations that can currently be counted are summarized below:

```@docs
Gap
Ambiguous
Certain
Match
Mismatch
Conserved
Mutated
Transition
Transversion
```

### The `count_sites` method

```@docs
count_sites
```

```julia
count_sites(Match, [dna"ATCGATCG", dna"ACCGATCG"])

count_sites(Mismatch, [dna"ATCGATCG", dna"ACCGATCG"])

count_sites(Conserved, [dna"ATCGATCG", dna"ACCGATCG"])

count_sites(Mutated, [dna"ATCGATCG", dna"ACCGATCG"])

count_sites(Transition, [rna"AUCGAUCG", rna"ACCGAUCG"])

count_sites(Transversion, [rna"AUCGAUCG", rna"ACCGAUCG"])
```


## File formats for representing genetic variation

This module supports some common file formats to read and write genetic
variations.


### VCF and BCF file formats

VCF and BCF file formats can be read using `VCFReader` and `BCFReader`,
respectively:
```julia
reader = open(VCFReader, "example.vcf")
for record in reader
    # do something
end
close(reader)
```

A reader first reads the header section of a file and creates a `VCFHeader`
object. The `header` function is provided to access the header object of a
reader:
```jlcon
julia> header(reader)
Bio.Var.VCFHeader:
  metainfo tags: fileformat fileDate source reference contig phasing INFO FILTER FORMAT
     sample IDs: NA00001 NA00002 NA00003

julia> find(header(reader), "FORMAT")
4-element Array{Bio.Var.VCFMetaInfo,1}:
 Bio.Var.VCFMetaInfo:
    tag: FORMAT
  value: ID="GT" Number="1" Type="String" Description="Genotype"
 Bio.Var.VCFMetaInfo:
    tag: FORMAT
  value: ID="GQ" Number="1" Type="Integer" Description="Genotype Quality"
 Bio.Var.VCFMetaInfo:
    tag: FORMAT
  value: ID="DP" Number="1" Type="Integer" Description="Read Depth"
 Bio.Var.VCFMetaInfo:
    tag: FORMAT
  value: ID="HQ" Number="2" Type="Integer" Description="Haplotype Quality"

```

`VCFMetaInfo` objects in the header support the following accessor functions:

| Accessor      | Description                          |
| :-------      | :----------                          |
| `metainfotag` | tag string                           |
| `metainfoval` | value string                         |
| `keys`        | keys of fields between '<' and '>'   |
| `values`      | values of fields between '<' and '>' |
| `[<key>]`     | value of a field with `key`          |

```jlcon
julia> metainfo = VCFMetaInfo("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
Bio.Var.VCFMetaInfo:
    tag: FORMAT
  value: ID="GT" Number="1" Type="String" Description="Genotype"

julia> metainfotag(metainfo)
"FORMAT"

julia> metainfoval(metainfo)
"<ID=GT,Number=1,Type=String,Description=\"Genotype\">"

julia> keys(metainfo)
4-element Array{String,1}:
 "ID"
 "Number"
 "Type"
 "Description"

julia> metainfo["ID"]
"GT"

```

`VCFRecord` and `BCFRecord` objects support the following accessor functions
(see the docstring of each accessor for the details):

| Accessor       | Description                    |
| :-------       | :----------                    |
| `chromosome`   | chromosome name                |
| `leftposition` | reference position             |
| `identifier`   | unique identifiers             |
| `reference`    | reference bases                |
| `alternate`    | alternate bases                |
| `quality`      | Phred-scaled quality score     |
| `filter_`      | filter status                  |
| `information`  | additional information         |
| `infokeys`     | keys of additional information |
| `format`       | genotype format                |
| `genotype`     | genotype information           |

```jlcon
julia> record = VCFRecord("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")

Bio.Var.VCFRecord:
   chromosome: 20
     position: 14370
   identifier: rs6054257
    reference: G
    alternate: A
      quality: 29.0
       filter: PASS
  information: NS=3 DP=14 AF=0.5 DB H2
       format: GT GQ DP HQ
     genotype: [1] 0|0 48 1 51,51 [2] 1|0 48 8 51,51

julia> chromosome(record)
Nullable{String}("20")

julia> leftposition(record)
Nullable{Int64}(14370)

julia> identifier(record)
1-element Array{String,1}:
 "rs6054257"

julia> reference(record)
Nullable{String}("G")

julia> alternate(record)
1-element Array{String,1}:
 "A"

julia> quality(record)
Nullable{Float64}(29.0)

julia> filter_(record)
1-element Array{String,1}:
 "PASS"

julia> information(record)
5-element Array{Pair{String,String},1}:
 "NS"=>"3"
 "DP"=>"14"
 "AF"=>"0.5"
 "DB"=>""
 "H2"=>""

julia> information(record, "AF")
"0.5"

julia> format(record)
4-element Array{String,1}:
 "GT"
 "GQ"
 "DP"
 "HQ"

julia> genotype(record)
2-element Array{Array{String,1},1}:
 String["0|0","48","1","51,51"]
 String["1|0","48","8","51,51"]

julia> genotype(record, 1)
4-element Array{String,1}:
 "0|0"
 "48"
 "1"
 "51,51"

julia> genotype(record, 1:2, "GT")
2-element Array{String,1}:
 "0|0"
 "1|0"

```

# MASH Distances

[MASH distances](http://doi.org/10.1186/s13059-016-0997-x), based on MinHash
sketches of genome sequences can provide rapid genome-scale sequence comparisons
when sequence distance (not specific mutations) are all that's required.

A MinHash sketch is made by taking the `s` smallest hash values for kmers of
length `k` for a given sequence. The genome distance for two genomes is then
essentially the [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index)
of the minhashes, with some additional modification to account for the size of
the kmers used.

You can generate a MinHash sketch using the `minhash()` function in `Bio.seq`.

```julia
using Bio.Seq

seq1 = dna"AAATAAGGCACAACTATGCAT"
sketch1 = minhash(seq, 5, 10)
```

Then, if you have MinHash sketches with the same parameters for two sequences,
you can determine the MASH distance between them.

```julia
seq2 = dna"AATTAACGCACGGACTGCGGTAAT"
sketch2 = minhash(seq, 5, 10)

using Bio.Var

mashdistance(sketch1, sketch2)
```

For more information on what size kmers and what size sketches are appropriate
for your use-case, see [Odnov et. al.](http://doi.org/10.1186/s13059-016-0997-x)
in _Genome Biology_. 
