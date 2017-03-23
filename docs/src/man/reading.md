# Reading and writing data

```@meta
CurrentModule = Bio
```

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
using Bio.Seq  # import FASTA
reader = open(FASTA.Reader, "hg38.fa")
# do something
close(reader)
```


## Reading by iteration

Readers in Bio.jl all read and return entries one at a time. The most convenient
way to do this by iteration:
```julia
reader = open(BED.Reader, "input.bed")
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
entry.  Some care is necessary when using this interface because `record` is
completely overwritten on each iteration:
```julia
reader = open(BED.Reader, "input.bed")
record = BED.Record()
while !eof(reader)
    read!(reader, record)
    # perform some operation on `record`
end
close(reader)
```

Empty record types that correspond to the file format be found using `eltype`,
making it easy to allocate an empty record for any reader stream:
```julia
record = eltype(stream)()
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
| GFF3        | `GFF3`   | `Bio.Intervals` | <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md> |
| bigBed      | `BigBed` | `Bio.Intervals` | <https://doi.org/10.1093/bioinformatics/btq351>                             |
| PDB         | `PDB`    | `Bio.Structure` | <http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html> |
| SAM         | `SAM`    | `Bio.Align`     | <https://samtools.github.io/hts-specs/SAMv1.pdf>                            |
| BAM         | `BAM`    | `Bio.Align`     | <https://samtools.github.io/hts-specs/SAMv1.pdf>                            |
| VCF         | `VCF`    | `Bio.Var`       | <https://samtools.github.io/hts-specs/VCFv4.3.pdf>                          |
| BCF         | `BCF`    | `Bio.Var`       | <https://samtools.github.io/hts-specs/VCFv4.3.pdf>                          |


### FASTA

* Reader type: `FASTA.Reader`
* Writer type: `FASTA.Writer`
* Element type: `FASTA.Record`

FASTA is a text-based file format for representing biological sequences. A FASTA
file stores a list of records with identifier, description, and sequence. The
template of a sequence record is:
```
>{identifier} {description}?
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
reader = open(FASTAReader, "sacCer.fa", index="sacCer.fa.fai")
chrIV = reader["chrIV"]  # directly read chromosome 4
```

```@docs
Bio.Seq.FASTA.Reader
Bio.Seq.FASTA.Writer
Bio.Seq.FASTA.Record
Bio.Seq.FASTA.identifier
Bio.Seq.FASTA.description
Bio.Seq.FASTA.sequence
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
# The default base quality encoding is Sanger.
reader = open(FASTQReader, "reads.fastq")
for record in reader
    # do something
end
close(reader)

# If performance is important, in-place reading will be much faster.
reader = open(FASTQReader, "reads.fastq")
record = FASTQSeqRecord{DNASequence}()
while !eof(reader)
    read!(reader, record)
    # do something
end
close(reader)
```

```@docs
Bio.Seq.FASTQReader
Bio.Seq.FASTQWriter
```


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
reader = open(TwoBitReader, "sacCer.2bit")  # load a random access index in the header
chrIV = reader["chrIV"]                     # directly read chromosome 4
```

```@docs
Bio.Seq.TwoBitReader
Bio.Seq.TwoBitWriter
```


### BED

* Reader type: `BED.Reader`
* Writer type: `BED.Writer`
* Element type: `BED.Record`

BED is a text-based file format for representing genomic annotations like genes,
transcripts, and so on. A BED file has tab-delimited and variable-length fields;
the first three fields denoting a genomic interval are mandatory.

This is an example of RNA transcripts:
```
chr9	68331023	68424451	NM_015110	0	+
chr9	68456943	68486659	NM_001206	0	-
```

```@docs
Bio.Intervals.BED.Reader
Bio.Intervals.BED.Writer
Bio.Intervals.BED.Record
Bio.Intervals.BED.chrom
Bio.Intervals.BED.chromstart
Bio.Intervals.BED.chromend
Bio.Intervals.BED.name
Bio.Intervals.BED.score
Bio.Intervals.BED.strand
Bio.Intervals.BED.thickstart
Bio.Intervals.BED.thickend
Bio.Intervals.BED.itemrgb
Bio.Intervals.BED.blockcount
Bio.Intervals.BED.blocksizes
Bio.Intervals.BED.blockstarts
```


### GFF3

* Reader type: `Bio.Intervals.GFF3.Reader`
* Writer type: `Bio.Intervals.GFF3.Writer`
* Element type: `Bio.Intervals.GFF3.Record`

GFF3 is a text-based file format for representing genomic annotations. The major
difference from BED is that is GFF3 is more structured and can include sequences
in the FASTA file format.

```@docs
Bio.Intervals.GFF3.Reader
Bio.Intervals.GFF3.Writer
```


### bigBed

* Reader type: `BigBedReader`
* Writer type: `BigBedWriter{T<:IO}`
* Element type: `Interval{BEDMetadata}` (alias: `BEDInterval`)

BigBed is a binary file format for representing genomic annotations and often
created from BED files. The bigBed files are indexed to quickly fetch specific
regions.

```@docs
Bio.Intervals.BigBedReader
Bio.Intervals.BigBedWriter
```


### PDB

PDB is a text-based file format for representing 3D macromolecular structures.
This has different reader interfaces from other file formats. Please consult the
[Bio.Structure](structure/) chapter for details.


### SAM

* Reader type: `SAMReader`
* Writer type: `SAMWriter{T<:IO}`
* Element type: `SAMRecord`

SAM is a text-based file format for representing sequence alignments.

```@docs
Bio.Align.SAMReader
Bio.Align.SAMWriter
```

### BAM

* Reader type: `BAMReader`
* Writer type: `BAMWriter`
* Element type: `BAMRecord`

BAM is a binary counterpart of the SAM file format.

When writing data in the BAM file format, the underlying output stream needs to
be wrapped with a `BGZFStream` object provided from
[BGZFStreams.jl](https://github.com/BioJulia/BGZFStreams.jl).

```@docs
Bio.Align.BAMReader
Bio.Align.BAMWriter
```


### VCF

* Reader type: `VCFReader`
* Writer type: `VCFWriter{T<:IO}`
* Element type: `VCFRecord`

VCF is a text-based file format for representing genetic variations.

```@docs
Bio.Var.VCFReader
Bio.Var.VCFWriter
```


### BCF

* Reader type: `BCFReader{T<:IO}`
* Writer type: `BCFWriter{T<:IO}`
* Element type: `BCFRecord`

BCF is a binary counterpart of the VCF file format.

```@docs
Bio.Var.BCFReader
Bio.Var.BCFWriter
```
