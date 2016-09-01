# Reading and writing data

Bio.jl has a unified interface for reading and writing files in a variety of
formats. Reader and writer type names have a prefix of the file format. For
example, files of a format `X` can be read using `XReader` and can be written
using `XWriter`.  To initialize a reader/writer of `X`, you can use one of the
following syntaxes:
```julia
# reader
open(::Type{XReader}, filepath::AbstractString, args...)
XReader(stream::IO, args...)

# writer
open(::Type{XWriter}, filepath::AbstractString, args...)
XWriter(stream::IO, args...)
```

For example, when reading a FASTA file, a reader for the FASTA file format can
be initialized as:
```julia
using Bio.Seq  # import FASTAReader
reader = open(FASTAReader, "hg38.fa")
# do something
close(reader)
```


## Reading by iteration

Readers in Bio.jl all read and return entries one at a time. The most convenient
way to do this by iteration:
```julia
reader = open(BEDReader, "input.bed")
for record in reader
    # perform some operation on entry
end
close(reader)
```


## In-place reading

Iterating through entries in a file is convenient, but for each entry in the
file, the reader must allocate, and ultimately the garbage collector must spend
time to deallocate it. For performance critical applications, a separate lower
level parsing interface can be used that avoid unnecessary allocation by
overwriting one entry. For files with a large number of small entries, this can
greatly speed up reading.

Instead of looping over a reader stream `read!` is called with a preallocated
entry.
```julia
reader = open(BEDReader, "input.bed")
record = BEDInterval()
while !eof(reader)
    read!(reader, record)
    # perform some operation on `entry`
end
close(reader)
```

Some care is necessary when using this interface. Because `entry` is completely
overwritten on each iteration, one must manually copy any field from `entry`
that should be preserved. For example, if we wish to save the `seqname` field
from `entry` when parsing BED, we must call `copy(entry.seqname)`.

Empty entry types that correspond to the file format be found using `eltype`,
making it easy to allocate an empty entry for any reader stream.

```julia
entry = eltype(stream)()
```


## Writing data

A FASTA file will be created as follows:
```julia
writer = open(FASTAWriter, "out.fa")
write(writer, FASTASeqRecord("seq1", dna"ACGTN"))
write(writer, FASTASeqRecord("seq2", dna"TTATA", "AT rich"))
close(writer)
```

Another way is using Julia's do-block syntax, which closes the data file after
finished writing:
```julia
open(FASTAWriter, "out.fa") do writer
    write(writer, FASTASeqRecord("seq1", dna"ACGTN"))
    write(writer, FASTASeqRecord("seq2", dna"TTATA", "AT rich"))
end
```


## Supported file formats

The following table summarizes supported file formats.

| File format | Prefix   | Module          | Specification                                                               |
| :---------- | :------- | :-----          | :------------                                                               |
| FASTA       | `FASTA`  | `Bio.Seq`       | <https://en.wikipedia.org/wiki/FASTA_format>                                |
| FASTQ       | `FASTQ`  | `Bio.Seq`       | <https://en.wikipedia.org/wiki/FASTQ_format>                                |
| .2bit       | `TwoBit` | `Bio.Seq`       | <http://genome.ucsc.edu/FAQ/FAQformat.html#format7>                         |
| BED         | `BED`    | `Bio.Intervals` | <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>                        |
| bigBed      | `BigBed` | `Bio.Intervals` | <https://doi.org/10.1093/bioinformatics/btq351>                             |
| PDB         | `PDB`    | `Bio.Structure` | <http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html> |


### FASTA

* Reader type: `FASTAReader{S<:Sequence}`
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
[fai](http://www.htslib.org/doc/faidx.html), the reader supports random access
to FASTA records, which would be useful when accessing specific parts of a huge
genome sequence:
```julia
reader = open(FASTAReader, "sacCer.fa")  # find and read "sacCer.fa.fai" file
chrIV = reader["chrIV"]  # directly read chromosome 4
```


### FASTQ

* Reader type: `FASTQReader{S<:Sequence}`
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

To read a file containing such records, one could use:
```julia
reader = open(FASTQReader, "reads.fastq", Seq.SANGER_QUAL_ENCODING)
seqrec = eltype(reader)()
while !eof(reader)
    read!(reader, seqrec)
    # ... process sequence
end
```

This assumes that the quality scores are in [Sanger
encoding](https://en.wikipedia.org/wiki/FASTQ_format#Encoding).


### .2bit

* Reader type: `TwoBitReader{T<:IO}`
* Writer type: `TwoBitWriter{T<:IO}`
* Element type: `SeqRecord{ReferenceSequence,Vector{UnitRange{Int}}}`

.2bit is a binary file format designed for storing a genome consists of multiple
chromosomal sequences. The reading speed is often an order of magnitude faster
than that of FASTA and the file size is smaller. However, since the .2bit file
format is specialized for genomic sequences, it cannot store either RNA or amino
acid sequences.

Like FASTA, the .2bit reader supports random access using an index included in
the header section of a .2bit file:
```julia
reader = open(TwoBitReader, "sacCer.2bit")  # read the header and load a random access index
chrIV = reader["chrIV"]  # directly read chromosome 4
```


### BED

* Reader type: `BEDReader`
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
This has different reader interfaces from other file formats. Please consult the
[Bio.Structure](structure/) chapter for details.
