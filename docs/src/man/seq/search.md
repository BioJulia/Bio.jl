```@meta
CurrentModule = Bio.Seq
DocTestSetup = quote
    using Bio.Seq
end
```

### Sequence search

Three kinds of on-line search functions are provided:

1. Exact search
2. Approximate search
3. Regular expression search

These are all specialized for biological sequences and ambiguities of symbols
are considered.

#### Exact search

Exact search functions search for an occurrence of the query symbol or
sequence. Four functions, `search`, `searchindex`, `rsearch`, and
`rsearchindex` are available:
```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> search(seq, DNA_G)  # search a query symbol
4:4

julia> query = dna"AGC";

julia> search(seq, query)  # search a query sequence
3:5

julia> searchindex(seq, query)
3

julia> rsearch(seq, query)  # similar to `search` but in the reverse direction
8:10

julia> rsearchindex(seq, query)  # similar to `searchindex` but in the reverse direction
8

```

These search functions take ambiguous symbols into account. That is, if two
symbols are compatible (e.g. `DNA_A` and `DNA_N`), they match when searching an
occurrence. In the following example, 'N' is a wild card that matches any
symbols:
```jldoctest
julia> search(dna"ACNT", DNA_N)  # 'A' matches 'N'
1:1

julia> search(dna"ACNT", dna"CGT")  # 'N' matches 'G'
2:4

julia> search(dna"ACGT", dna"CNT")  # 'G' matches 'N'
2:4

```

The exact sequence search needs preprocessing phase of query sequence before
searching phase. This would be enough fast for most search applications. But
when searching a query sequence to large amounts of target sequences, caching
the result of preprocessing may save time. The `ExactSearchQuery` creates such
a preprocessed query object and is applicable to the search functions:
```jldoctest
julia> query = ExactSearchQuery(dna"ATT");

julia> search(dna"ATTTATT", query)
1:3

julia> rsearch(dna"ATTTATT", query)
5:7

```


#### Approximate search

The approximate search is similar to the exact search but allows a specific
number of errors. That is, it tries to find a subsequence of the target sequence
within a specific [Levenshtein
distance](https://en.wikipedia.org/wiki/Levenshtein_distance) of the query
sequence:
```jldoctest
julia> seq = dna"ACAGCGTAGCT";

julia> approxsearch(seq, dna"AGGG", 0)  # nothing matches with no errors
0:-1

julia> approxsearch(seq, dna"AGGG", 1)  # seq[3:5] matches with one error
3:6

julia> approxsearch(seq, dna"AGGG", 2)  # seq[1:4] matches with two errors
1:4

```

Like the exact search functions, four kinds of functions (`approxsearch`,
`approxsearchindex`, `approxrsearch`, and `approxrsearchindex`) are available:
```jldoctest
julia> seq = dna"ACAGCGTAGCT"; pat = dna"AGGG";

julia> approxsearch(seq, pat, 2)        # return the range (forward)
1:4

julia> approxsearchindex(seq, pat, 2)   # return the starting index (forward)
1

julia> approxrsearch(seq, pat, 2)       # return the range (backward)
8:11

julia> approxrsearchindex(seq, pat, 2)  # return the starting index (backward)
8

```

Preprocessing can be cached in an `ApproximateSearchQuery` object:
```jldoctest
julia> query = ApproximateSearchQuery(dna"AGGG");

julia> approxsearch(dna"AAGAGG", query, 1)
2:5

julia> approxsearch(dna"ACTACGT", query, 2)
4:6

```

### Regular expression search

Query patterns can be described in regular expressions. The syntax supports
a subset of Perl and PROSITE's notation.

The Perl-like syntax starts with `biore` (**bio**logical **re**gular expression)
and ends with a symbol option: "dna", "rna" or "aa". For example, `biore"A+"dna`
is a regular expression for DNA sequences and `biore"A+"aa` is for amino acid
sequences. The symbol options can be abbreviated to its first character: "d",
"r" or "a", respectively.

Here are examples of using the regular expression for `BioSequence`s:
```jldoctest
julia> match(biore"A+C*"dna, dna"AAAACC")
Nullable{Bio.Seq.RE.RegexMatch{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}}}}(RegexMatch("AAAACC"))

julia> match(biore"A+C*"d, dna"AAAACC")
Nullable{Bio.Seq.RE.RegexMatch{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}}}}(RegexMatch("AAAACC"))

julia> ismatch(biore"A+C*"dna, dna"AAC")
true

julia> ismatch(biore"A+C*"dna, dna"C")
false

```

`match` always returns a `Nullable` object and it should be null if no match is
found.

The table below summarizes available syntax elements.

| Syntax | Description | Example |
|:------:|:------------|:--------|
| `\|` | alternation | `"A\|T"` matches `"A"` and `"T"` |
| `*` | zero or more times repeat | `"TA*"` matches `"T"`, `"TA"` and `"TAA"` |
| `+` | one or more times repeat | `"TA+"` matches `"TA"` and `"TAA"` |
| `?` | zero or one time | `"TA?"` matches `"T"` and `"TA"` |
| `{n,}` | `n` or more times repeat | `"A{3,}"` matches `"AAA"` and `"AAAA"` |
| `{n,m}` | `n`-`m` times repeat | `"A{3,5}"` matches `"AAA"`, `"AAAA"` and `"AAAAA"`|
| `^` | the start of the sequence | `"^TAN*"` matches `"TATGT"` |
| `$` | the end of the sequence | `"N*TA$"` matches `"GCTA"` |
| `(...)` | pattern grouping | `"(TA)+"` matches `"TA"` and `"TATA"` |
| `[...]` | one of symbols | `"[ACG]+"` matches `"AGGC"` |

`eachmatch`, `matchall`, and `search` are also defined like usual strings:
```jldoctest
julia> matchall(biore"TATA*?"d, dna"TATTATAATTA")  # overlap (default)
4-element Array{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},1}:
 TAT  
 TAT  
 TATA
 TATAA

julia> matchall(biore"TATA*"d, dna"TATTATAATTA", false)  # no overlap
2-element Array{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},1}:
 TAT  
 TATAA

julia> search(dna"TATTATAATTA", biore"TATA*"d)
1:3

julia> search(dna"TATTATAATTA", biore"TATA*"d, 2)
4:8

```

Notewothy differences from strings are:

* Ambiguous characters match any compatible characters (e.g. `biore"N"d` is equivalent to `biore"[ACGT]"d`).
* Whitespaces are ignored (e.g. `biore"A C G"d` is equivalent to `biore"ACG"d`).

The PROSITE notation is described in [ScanProsite - user
manual](http://prosite.expasy.org/scanprosite/scanprosite_doc.html). The syntax
supports almost all notations including the extended syntax. The PROSITE
notation starts with `prosite` prefix and no symbol option is needed because it
always describes patterns of amino acid sequences:
```jldoctest
julia> match(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVPQARG")
Nullable{Bio.Seq.RE.RegexMatch{Bio.Seq.BioSequence{Bio.Seq.AminoAcidAlphabet}}}(RegexMatch("CPVPQARG"))

julia> match(prosite"[AC]xVx(4){ED}", aa"CPVPQARG")
Nullable{Bio.Seq.RE.RegexMatch{Bio.Seq.BioSequence{Bio.Seq.AminoAcidAlphabet}}}(RegexMatch("CPVPQARG"))

```


## Sequence composition

Sequence composition can be easily calculated using the `composition` function:
```jldoctest
julia> comp = composition(dna"ACGAG")
DNA Composition:
  DNA_Gap => 0
  DNA_A   => 2
  DNA_C   => 1
  DNA_M   => 0
  DNA_G   => 2
  DNA_R   => 0
  DNA_S   => 0
  DNA_V   => 0
  DNA_T   => 0
  DNA_W   => 0
  DNA_Y   => 0
  DNA_H   => 0
  DNA_K   => 0
  DNA_D   => 0
  DNA_B   => 0
  DNA_N   => 0

julia> comp[DNA_A]
2

julia> comp[DNA_T]
0

```

To accumulate composition statistics of multiple sequences, `merge!` can be used
as follows:
```@repl
# initiaize an empty composition counter
comp = composition(dna"");

# iterate over sequences and accumulate composition statistics into `comp`
for seq in seqs
    merge!(comp, composition(seq))
end

# or functional programming style in one line
foldl((x, y) -> merge(x, composition(y)), composition(dna""), seqs)
```

`composition` is also applicable to a *k*-mer iterator:
```jldoctest
julia> comp = composition(each(DNAKmer{4}, dna"ACGT"^100));

julia> comp[DNAKmer("ACGT")]
100

julia> comp[DNAKmer("CGTA")]
99

```
