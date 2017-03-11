## Sequence records

The `SeqRecord` type is used to represent a named sequence, optionally with
accompanying metadata. It is defined as:
sequence. It is as:
```jlcon
type SeqRecord{S, T}
    name::String
    seq::S
    metadata::T
end
```

The type of the `metadata` field depends on the source of the sequence record.
For example, if a record is read from a FASTA file, metadata contains the
description field. If from a FASTQ file, a quality scores assigned to base calls
during sequencing.

The following accessors are defined for the `SeqRecord` type:
```@docs
seqname
sequence
metadata
```
