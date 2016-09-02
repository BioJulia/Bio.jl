# Intervals: Genomic Interval Manipulation

```@meta
CurrentModule = Bio.Intervals
```

The `Intervals` module consists of tools for working efficiently with genomic
intervals.


## Interval types

Intervals in Bio.jl are consistent with ranges in Julia: *1-based and
end-inclusive*. When data is read from formats with different representations
(i.e. 0-based and/or end-exclusive) they are always converted automatically.
Similarly when writing data. You should not have to reason about off-by-one
errors due to format differences while using functionality provided in Bio.jl.

The `Interval` type is defined as
```julia
type Interval{T} <: AbstractInterval{Int64}
    seqname::StringField
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end
```

Similarly to the `SeqRecord` type in the `Seq` module, `Interval` is
parameterized on metadata type, which lets it efficiently and precisely
be specialized to represent intervals from a variety of formats.

Strand is represented by the `Strand` type which can take four possible values:
```julia
STRAND_NA   # strand is unknown or inapplicable
STRAND_POS  # positive strand
STRAND_NEG  # negative strand
STRAND_BOTH # non-strand-specific feature
```

## Collections of intervals

Collections of intervals are represented using the `IntervalCollection` type,
which is a general purpose indexed container for intervals. It supports fast
intersection operations as well as insertion, deletion, and sorted iteration.

Interval collections can be initialized by inserting elements one by one using
`push!`.

```julia
# The type parameter (Void here) indicates the interval metadata type.
incol = IntervalCollection{Void}()

for i in 1:100:10000
    push!(incol, Interval("chr1", i, i + 99, STRAND_POS, nothing))
end
```

Incrementally building an interval collection like this works, but
`IntervalCollection` also has a bulk insertion constructor that is able to build
the indexed data structure extremely efficiently from an array of intervals.

```julia
incol = IntervalCollection([Interval("chr1", i, i + 99, STRAND_POS, nothing) for i in 1:100:10000])
```

Bulding `IntervalCollections` in one shot like this should be preferred when
it's convenient or speed in an issue.


## Intersection

There are number of `intersect` function in the Intervals module. They follow
two patterns: interval versus collection queries which return an iterator over
intervals in the collection that intersect the query, and collection versus
collection queries which iterate over all pairs of intersecting intervals.

```@docs
intersect
```


## Interval streams

Intervals need not necessarily stored in an indexed data structure for efficient
intersection to be practical. Two collections of intervals need only be both
sorted to compute all intersecting pairs. This is particularly useful in
genomics where datasets are sometimes so large that loading them entirely into
memory is not practical.

The Intervals module is able to intersect any two iterators that yield intervals
in sorted order, which we refer to as "interval streams". An
`IntervalCollection` is also a interval stream, but so is a sorted array of
intervals, and parsers over interval file formats. This allows for a very
general notion of intersection.

```julia
for (x, y) in intersect(open(BEDReader, "x_features.bed"), open(BEDReader, "y_features.bed"))
    println("Intersection found between ", x, " and ", y)
end
```

An exception will be thrown if an interval in encountered out of order while
processing an interval stream. Ordering of intervals has one complication: there
is not necessarily a standardized way to order sequence names. By default in
Bio.jl intervals are sorted using a special `alphanum_isless` comparison
function that compares numbers numerically if they exist in string, so that
names like `chr1, chr2, chr10` end up in their natural order.

The `intersect` function takes as an optional parameter an `isless` function to
use to compare sequence names to account for arbitrary sequence name orderings.

```julia
# assume lexigraphic ordering for sequence names
for (x, y) in intersect(xs, ys, isless)
    println("Intersection found between ", a, " and ", b)
end
```

A special sort of intersection can also be performed on a `IntervalStreams`
against itself to produce "coverage intervals".

```@docs
coverage
```


## Random access of intervals

It is often the case that the user is interested in specific regions of a file.
In such a case, reading all data from top to the specific regions is a
time-wasting process when handing a big data file.

[Tabix](http://www.htslib.org/doc/tabix.html) is a tool to index an interval
data file and makes it possible to quickly retrieve data overlapping with a
specific genomic region. Indexing a file using the tabix tool requires that the
data file is sorted and then compressed in a special file format called
[BGZF](https://github.com/BioJulia/BGZFStreams.jl). The BGZF file format is a
gzip-compliant data format and readable using standard command-line tools like
`gzip`. However, creating a BGZF file needs a purpose-built tool called `bgzip`.
Once a BGZF file and its index file are prepared using `bgzip` and `tabix`, the
user can iterate over intervals overlapping with a specified region:
```julia
reader = open(BEDReader, "data.bed.gz")
for interval in intersect(reader, "chr2", 1_000_000:2_000_000)
    # do something
end
```

Internally, this random accessing function uses an index type named
`Tabix`.  If the user knows well about the BGZF file format, it is also possible
to use its interfaces to read a formatted file that are not supported in Bio.jl:
```julia
index = Tabix("data.xxx.gz.tbi")
stream = BGZFStream("data.xxx.gz")
for chunk in overlapchunks(index, "chr2", 1_000_000:2_000_000)
    seekstart(stream, chunk.start)
    for line in eachline(stream)
        values = split(chomp(line), '\t')
        seqname = values[index.columns[1]]
        start = parse(Int, values[index.columns[2]])
        stop = parse(Int, values[index.columns[3]])
        if seqname == "chr2" && intersect(start:stop, 1_000_000:2_000_000)
            # do something
        end
    end
end
```
