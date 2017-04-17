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
writer = open(FASTA.Writer, "out.fa")
write(writer, FASTA.Record("seq1", dna"ACGTN"))
write(writer, FASTA.Record("seq2", "AT rich", dna"TTATA"))
close(writer)
```

Another way is using Julia's do-block syntax, which closes the data file after
finished writing:
```julia
open(FASTA.Writer, "out.fa") do writer
    write(writer, FASTA.Record("seq1", dna"ACGTN"))
    write(writer, FASTA.Record("seq2", "AT rich", dna"TTATA"))
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

* Reader type: `FASTQ.Reader`
* Writer type: `FASTQ.Writer`
* Element type: `FASTQ.Record`

FASTQ is a text-based file format for representing DNA sequences along with
qualities for each base. A FASTQ file stores a list of sequence records in the
following format:
```
@{identifier} {description}?
{sequence}
+
{quality}
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
reader = open(FASTQ.Reader, "reads.fastq")
for record in reader
    # do something
end
close(reader)

# If performance is important, in-place reading will be much faster.
reader = open(FASTQ.Reader, "reads.fastq")
record = FASTQ.Record()
while !eof(reader)
    read!(reader, record)
    # do something
end
close(reader)
```

```@docs
Bio.Seq.FASTQ.Reader
Bio.Seq.FASTQ.Writer
Bio.Seq.FASTQ.Record
Bio.Seq.FASTQ.identifier
Bio.Seq.FASTQ.description
Bio.Seq.FASTQ.sequence
Bio.Seq.FASTQ.quality
```


### .2bit

* Reader type: `TwoBit.Reader{T<:IO}`
* Writer type: `TwoBit.Writer{T<:IO}`
* Element type: `TwoBit.Record`

.2bit is a binary file format designed for storing a genome consists of multiple
chromosomal sequences. The reading speed is often an order of magnitude faster
than that of FASTA and the file size is smaller. However, since the .2bit file
format is specialized for genomic sequences, it cannot store either RNA or amino
acid sequences.

Like FASTA, the .2bit reader supports random access using an index included in
the header section of a .2bit file:
```julia
reader = open(TwoBit.Reader, "sacCer.2bit")  # load a random access index in the header
chrIV = reader["chrIV"]                      # directly read chromosome 4
```

```@docs
Bio.Seq.TwoBit.Reader
Bio.Seq.TwoBit.Writer
Bio.Seq.TwoBit.Record
Bio.Seq.TwoBit.seqnames
Bio.Seq.TwoBit.sequence
Bio.Seq.TwoBit.maskedblocks
```

### ABIF

* Reader type: `AbifReader{T<:IO}`
ABIF is a binary file format for storing data produced by sequencers, those developed by Applied Biosystems, Inc.
When the file is opened, we save all the existing tags, so we can read only the tags that are needed.
```julia
reader = open(AbifReader, "3100.ab1")       # load a random access
data   = reader["DATA"]                     # directly read all existing `DATA` Tags
data   = reader[1]                          # directly read Tag at index
data  = tags(reader)                       # return all existing tags

# iterator by all tags
for (key, value) in getindex(reader, data)
end
```

```@docs
Bio.Seq.AbifReader
Bio.Seq.tags
Bio.Seq.elements
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
Bio.Intervals.GFF3.Record
Bio.Intervals.GFF3.isfeature
Bio.Intervals.GFF3.isdirective
Bio.Intervals.GFF3.iscomment
Bio.Intervals.GFF3.seqid
Bio.Intervals.GFF3.source
Bio.Intervals.GFF3.featuretype
Bio.Intervals.GFF3.seqstart
Bio.Intervals.GFF3.seqend
Bio.Intervals.GFF3.score
Bio.Intervals.GFF3.strand
Bio.Intervals.GFF3.phase
Bio.Intervals.GFF3.attributes
Bio.Intervals.GFF3.content
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

* Reader type: `SAM.Reader`
* Writer type: `SAM.Writer`
* Element type: `SAM.Record`

SAM is a text-based file format for representing sequence alignments.

```@docs
Bio.Align.SAM.Reader
Bio.Align.SAM.header

Bio.Align.SAM.Header
Base.find(header::Bio.Align.SAM.Header, key::AbstractString)

Bio.Align.SAM.Writer

Bio.Align.SAM.MetaInfo
Bio.Align.SAM.iscomment
Bio.Align.SAM.tag
Bio.Align.SAM.value
Bio.Align.SAM.keyvalues

Bio.Align.SAM.Record
Bio.Align.SAM.flag
Bio.Align.SAM.ismapped
Bio.Align.SAM.refname
Bio.Align.SAM.position
Bio.Align.SAM.rightposition
Bio.Align.SAM.isnextmapped
Bio.Align.SAM.nextrefname
Bio.Align.SAM.nextposition
Bio.Align.SAM.mappingquality
Bio.Align.SAM.cigar
Bio.Align.SAM.alignment
Bio.Align.SAM.alignlength
Bio.Align.SAM.tempname
Bio.Align.SAM.templength
Bio.Align.SAM.sequence
Bio.Align.SAM.seqlength
Bio.Align.SAM.quality
Bio.Align.SAM.auxdata
```

This module provides 16-bit flags defined in the SAM specs:

| Flag                      | Bit       | Description                                                        |
| :------------------------ | :-------- | :----------------------------------------------------------------- |
| `SAM.FLAG_PAIRED`         | `0x0001`  | template having multiple segments in sequencing                    |
| `SAM.FLAG_PROPER_PAIR`    | `0x0002`  | each segment properly aligned according to the aligner             |
| `SAM.FLAG_UNMAP`          | `0x0004`  | segment unmapped                                                   |
| `SAM.FLAG_MUNMAP`         | `0x0008`  | next segment in the template unmapped                              |
| `SAM.FLAG_REVERSE`        | `0x0010`  | SEQ being reverse complemented                                     |
| `SAM.FLAG_MREVERSE`       | `0x0020`  | SEQ of the next segment in the template being reverse complemented |
| `SAM.FLAG_READ1`          | `0x0040`  | the first segment in the template                                  |
| `SAM.FLAG_READ2`          | `0x0080`  | the last segment in the template                                   |
| `SAM.FLAG_SECONDARY`      | `0x0100`  | secondary alignment                                                |
| `SAM.FLAG_QCFAIL`         | `0x0200`  | not passing filters, such as platform/vendor quality controls      |
| `SAM.FLAG_DUP`            | `0x0400`  | PCR or optical duplicate                                           |
| `SAM.FLAG_SUPPLEMENTARY`  | `0x0800`  | supplementary alignment                                            |


### BAM

* Reader type: `BAM.Reader`
* Writer type: `BAM.Writer`
* Element type: `BAM.Record`

BAM is a binary counterpart of the [SAM](@ref) file format.

When writing data in the BAM file format, the underlying output stream needs to
be wrapped with a `BGZFStream` object provided from
[BGZFStreams.jl](https://github.com/BioJulia/BGZFStreams.jl).

Flags and the header type are defined in the [SAM](@ref) module.

```@docs
Bio.Align.BAM.Reader
Bio.Align.BAM.header

Bio.Align.BAM.Writer

Bio.Align.BAM.Record
Bio.Align.BAM.flag
Bio.Align.BAM.ismapped
Bio.Align.BAM.refid
Bio.Align.BAM.refname
Bio.Align.BAM.position
Bio.Align.BAM.rightposition
Bio.Align.BAM.isnextmapped
Bio.Align.BAM.nextrefid
Bio.Align.BAM.nextrefname
Bio.Align.BAM.nextposition
Bio.Align.BAM.mappingquality
Bio.Align.BAM.cigar
Bio.Align.BAM.cigar_rle
Bio.Align.BAM.alignment
Bio.Align.BAM.alignlength
Bio.Align.BAM.tempname
Bio.Align.BAM.templength
Bio.Align.BAM.sequence
Bio.Align.BAM.seqlength
Bio.Align.BAM.quality
Bio.Align.BAM.auxdata
```


### VCF

* Reader type: `VCF.Reader`
* Writer type: `VCF.Writer{T<:IO}`
* Element type: `VCF.Record`

VCF is a text-based file format for representing genetic variations.

```@docs
Bio.Var.VCF.Reader
Bio.Var.VCF.Writer
Bio.Var.VCF.Record
Bio.Var.VCF.chrom
Bio.Var.VCF.pos
Bio.Var.VCF.id
Bio.Var.VCF.ref
Bio.Var.VCF.alt
Bio.Var.VCF.qual
Bio.Var.VCF.filter
Bio.Var.VCF.info
Bio.Var.VCF.format
```


### BCF

* Reader type: `BCF.Reader{T<:IO}`
* Writer type: `BCF.Writer{T<:IO}`
* Element type: `BCF.Record`

BCF is a binary counterpart of the VCF file format.

```@docs
Bio.Var.BCF.Reader
Bio.Var.BCF.Writer
Bio.Var.BCF.Record
Bio.Var.BCF.chrom
Bio.Var.BCF.pos
Bio.Var.BCF.rlen
Bio.Var.BCF.qual
Bio.Var.BCF.ref
Bio.Var.BCF.alt
Bio.Var.BCF.filter
Bio.Var.BCF.info
Bio.Var.BCF.genotype
```
