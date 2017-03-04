# Bio.Var: Biological Variation.

```@meta
CurrentModule = Bio.Var
DocTestSetup = quote
    using Bio.Seq
    using Bio.Var
end
```

## Identifying and counting mutations

You can count the numbers of different types of mutations in a pairwise manner
for a set of nucleotide sequences.

### Different types of mutation

The types of mutations that can currently be counted are `AnyMutation`s,
`TransitionMutation`s, and `TransversionMutation`s.

```@docs
AnyMutation
TransitionMutation
TransversionMutation
```

### The `count_mutations` method

Mutations are counted using the `count_mutations` method.
The method outputs a tuple. The first value is the number of mutations counted.
The second value is the number of sites examined. Sites which have gaps and
uncertain nucleotides are not examined and so this second value will be less
than the length of the two biological sequences.

```jldoctest
julia> count_mutations(AnyMutation, [dna"ATCGATCG", dna"ACCGATCG"])
([1],[8])

julia> count_mutations(TransitionMutation, [dna"ATCGATCG", dna"ACCGATCG"])
([1],[8])

julia> count_mutations(TransversionMutation, [dna"ATCGATCG", dna"ACCGATCG"])
([0],[8])

julia> count_mutations(TransitionMutation, TransversionMutation, [dna"ATCGATCG", dna"ACCGATCG"])
([1],[0],[8])

julia> count_mutations(AnyMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])
([1],[8])

julia> count_mutations(TransitionMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])
([1],[8])

julia> count_mutations(TransversionMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])
([0],[8])

julia> count_mutations(TransitionMutation, TransversionMutation, [rna"AUCGAUCG", rna"ACCGAUCG"])
([1],[0],[8])
```

### The `is_mutation` method

```@docs
is_mutation
```

## Computing evolutionary and genetic distances

Just as you can count the number of mutations between two nucleotide sequences,
you can compute the evolutionary distance between two nucleotide sequences.

### Different evolutionary distance measures

The types of distances that can currently be computed are described below.

```@docs
Count{T}
Proportion{T}
JukesCantor69
Kimura80
```

### The distance method

```@docs
distance
distance{T<:MutationType,A<:NucleicAcidAlphabet}(::Type{Count{T}}, seqs::Vector{BioSequence{A}})
distance{T<:MutationType,A<:NucleicAcidAlphabet}(::Type{Proportion{T}}, seqs::Vector{BioSequence{A}})
distance{A<:NucleicAcidAlphabet}(::Type{JukesCantor69}, seqs::Vector{BioSequence{A}})
distance{A<:NucleicAcidAlphabet}(::Type{Kimura80}, seqs::Vector{BioSequence{A}})
```


## File formats for represeting genetic variations

The module suports some common file formats to read and write genetic
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
(see docstring of each accessor for the details):

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
"G"

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

```
