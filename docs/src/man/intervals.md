# Intervals: Genomic Interval Manipulation

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
    seqname::ASCIIString
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end
```

Similarly to the `ReqRecord` type in the `Seq` module, `Interval` is
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
# The type parameter (Nothing here) indicates the interval metadata type.
incol = IntervalCollection{Nothing}()

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

    {docs}
    intersect


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
for (x, y) in intersect(open("x_features.bed", BED), open("y_features.bed", BED))
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

    {docs}
    coverage

