# Reading and writing data

Bio.jl has a unified interface for reading and writing files in a variety of
formats. To initialize a parser for a particular format, the `open` method is
extended with a file format type parameter (currently supported formats types
are `BED`, `FASTQ`, and `FASTA`).

```julia
open{T <: FileFormat}(filename::String, ::Type{T}; memory_map::Bool=false)
open{T <: FileFormat}(source::IO, ::Type{T})
open{T <: FileFormat}(data::Vector{UInt8}, ::Type{T})
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

