# Reading and writing data

Bio.jl has a unified interface for reading and writing files in a variety of
formats. To initialize a parser for a particular format, the `open` method is
extended with a file format type parameter `T`:
```julia
open{T<:FileFormat}(filename::AbstractString, ::Type{T})
open{T<:FileFormat}(source::IO, ::Type{T})
open{T<:FileFormat}(data::Vector{UInt8}, ::Type{T})
```

For example, when reading a FASTA file, a parser for the FASTA file format can
be initialized as:
```julia
using Bio.Seq  # import FASTA
parser = open("hg38.fa", FASTA)
```


## Parsing by iteration

Parsers in Bio.jl all read and return entries one at a time. The most convenient
way to do this by iteration.

```julia
stream = open("input.bed", BED)
for entry in stream
    # perform some operation on entry
end
```


## In-place parsing

Iterating through entries in a file is convenient, but for each entry in the
file, the parser must allocate, and ultimately the garbage collector must spend
time to deallocate it. For performance critical applications, a separate lower
level parsing interface can be used that avoid unnecessary allocation by
overwriting one entry. For files with a large number of small entries, this can
greatly speed up reading.

Instead of looping over a parser stream `read!` is called with a preallocated
entry.
```julia
stream = open("input.bed", BED)
entry = BEDInterval()
while !eof(stream)
    read!(input, entry)
    # perform some operation on `entry`
end
```

Some care is necessary when using this interface. Because `entry` is completely
overwritten on each iteration, one must manually copy any field from `entry`
that should be preserved. For example, if we wish to save the `seqname` field
from `entry` when parsing BED, we must call `copy(entry.seqname)`.

Empty entry types that correspond to the file format be found using `eltype`,
making it easy to allocate an empty entry for any parser stream.

```julia
entry = eltype(stream)()
```


## Writing data

Writing data into a stream has a uniform interface consistent with parsers. The
following code is a template of formatted serialization into a file:
```julia
# open a file of a particular file format in writing mode
out = open(<filepath>, "w", <format>)
# write a record into it
write(out, <record>)
# finally close it
close(out)
```

For example, a FASTA file will be created as follows:
```julia
out = open("out.fasta", "w", FASTA)
write(out, FASTASeqRecord("seq1", dna"ACGTN"))
write(out, FASTASeqRecord("seq2", dna"TTATATTATTGTAAA", "AT rich"))
# and more records
close(out)
```


## Supported file formats

The following table summarizes supported file formats.

| File format | Type parameter | Module          | Specification                                                               |
| :---------- | :------------- | :-----          | :------------                                                               |
| FASTA       | `FASTA`        | `Bio.Seq`       | <https://en.wikipedia.org/wiki/FASTA_format>                                |
| FASTQ       | `FASTQ`        | `Bio.Seq`       | <https://en.wikipedia.org/wiki/FASTQ_format>                                |
| .2bit       | `TwoBit`       | `Bio.Seq`       | <http://genome.ucsc.edu/FAQ/FAQformat.html#format7>                         |
| BED         | `BED`          | `Bio.Intervals` | <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>                        |
| bigBed      | `BigBed`       | `Bio.Intervals` | <https://doi.org/10.1093/bioinformatics/btq351>                             |
| PDB         | `PDB`          | `Bio.Structure` | <http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html> |


### FASTA

* Parser type: `FASTAParser{S<:Sequence}`
* Writer type: `FASTAWriter{T<:IO}`
* Element type: `SeqRecord{S,FASTAMetadata}` (alias: `FASTASeqRecord{S}`)

FASTA is a text-based file format for representing biological sequences. A
FASTA file stores a list of sequence records with name, description, and
sequence. The template of a sequence record is:
```
>{name} {description}?
{sequence}
```

Here is an example of a chromosomal sequence:
```
>chrI chromosome 1
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACC
CACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG
```

Usually sequence records will be read sequentially from a file by iteration.
But if the FASTA file has an auxiliary index file formatted in
[fai](http://www.htslib.org/doc/faidx.html), the parser supports random access
to FASTA records, which would be useful when accessing specific parts of a huge
genome sequence:
```julia
parser = open("sacCer.fa", FASTA)  # find and read "sacCer.fa.fai" file
chrIV = parser["chrIV"]  # directly read chromosome 4
```


### FASTQ

* Parser type: `FASTQParser{S<:Sequence}`
* Writer type: `FASTQWriter{T<:IO}`
* Element type: `SeqRecord{S,FASTQMetadata}` (alias: `FASTQSeqRecord{S}`)

FASTQ is a text-based file format for representing DNA sequences along with
qualities for each base. A FASTQ file stores a list of sequence records in the
following format:
```
@{name} {description}?
{sequence}
+
{qualities}
```

Here is an example of one record from a FASTQ file:
```
@FSRRS4401BE7HA
tcagTTAAGATGGGAT
+
###EEEEEEEEE##E#
```

To parse a file containing such records, one could use:
```julia
parser = open("reads.fastq", FASTQ, Seq.SANGER_QUAL_ENCODING)
seqrec = eltype(parser)()
while !eof(parser)
    read!(parser, seqrec)
    # ... process sequence
end
```

This assumes that the quality scores are in [Sanger
encoding](https://en.wikipedia.org/wiki/FASTQ_format#Encoding).


### .2bit

* Parser type: `TwoBitParser{T<:IO}`
* Writer type: `TwoBitWriter{T<:IO}`
* Element type: `SeqRecord{ReferenceSequence,Vector{UnitRange{Int}}}`

.2bit is a binary file format designed for storing a genome consists of multiple
chromosomal sequences. The reading speed is often an order of magnitude faster
than that of FASTA and the file size is smaller. However, since the .2bit file
format is specialized for genomic sequences, it cannot store either RNA or amino
acid sequences.

Like FASTA, the .2bit parser supports random access using an index included in
the header section of a .2bit file:
```julia
parser = open("sacCer.2bit", TwoBit)  # parse the header and load a random access index
chrIV = parser["chrIV"]  # directly read chromosome 4
```


### BED

* Parser type: `BEDParser`
* Writer type: `BEDWriter{T<:IO}`
* Element type: `Interval{BEDMetadata}` (alias: `BEDInterval`)

BED is a text-based file format for representing genomic annotations like genes,
transcripts, and so on. A BED file has tab-delimited and variable-length fields;
the first three fields denoting a genomic interval are mandatory.

This is an example of RNA transcripts:
```
chr9	68331023	68424451	NM_015110	0	+
chr9	68456943	68486659	NM_001206	0	-
```


### bigBed

BigBed is a binary file format for representing genomic annotations and often
created from BED files. The bigBed files are indexed to quickly fetch specific
regions.


### PDB

PDB is a text-based file format for representing 3D macromolecular structures.
This has different parser interfaces from other file formats. Please consult the
[Bio.Structure](structure/) chapter for details.
