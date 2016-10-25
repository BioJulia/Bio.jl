# Bio.Services: APIs for Web Services

```@meta
CurrentModule = Bio.Services
```

The `Bio.Services` module provides APIs for various web services.


## E-Utilities

E-Utilities provide a interface to Entrez databases at
[NCBI](https://www.ncbi.nlm.nih.gov/).  The APIs are defined in the
`Bio.Services.EUtils` module, which exports nine functions to access its databases:

| Function    | Description                                                                |
| :-------    | :----------                                                                |
| `einfo`     | Retrieve a list of databases or statistics for a database.                 |
| `esearch`   | Retrieve a list of UIDs matching a text query.                             |
| `epost`     | Upload or append a list of UIDs to the Entrez History server.              |
| `esummary`  | Retrieve document summaries for a list of UIDs.                            |
| `efetch`    | Retrieve formatted data records for a list of UIDs.                        |
| `elink`     | Retrieve UIDs linked to an input set of UIDs.                              |
| `egquery`   | Retrieve the number of available records in all databases by a text query. |
| `espell`    | Retrieve spelling suggestions.                                             |
| `ecitmatch` | Retrieve PubMed IDs that correspond to a set of input citation strings.    |

["The Nine E-utilities in
Brief"](https://www.ncbi.nlm.nih.gov/books/NBK25497/#_chapter2_The_Nine_Eutilities_in_Brief_)
summarizes all of the server-side programs corresponding to each function.

In this package, queries for databases are controlled by keyword parameters. For
example, some functions take `db` parameter to specify the target database.
Functions listed above take these parameters as keyword arguments and return
a `Response` object as follows:
```jlcon
julia> using Bio.Services.EUtils      # import the nine functions above

julia> res = einfo(db="pubmed")       # retrieve statistics of the PubMed database
Response(200 OK, 18 headers, 27360 bytes in body)

julia> write("pubmed.xml", res.data)  # save retrieved data into a file
27360

shell> head result.xml
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD einfo 20130322//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20130322/einfo.dtd">
<eInfoResult>
        <DbInfo>
        <DbName>pubmed</DbName>
        <MenuName>PubMed</MenuName>
        <Description>PubMed bibliographic record</Description>
        <DbBuild>Build161024-2207m.1</DbBuild>
        <Count>26590895</Count>
        <LastUpdate>2016/10/25 02:06</LastUpdate>

```

Let's see a few more examples of parameters.  The `term` parameter specifies a
search string (e.g. `esearch(db="gene", term="tumor AND human[ORGN]")`).  The
`id` parameter specifies a UID (or accession number) or a list of UIDs (e.g.
`efetch(db="protein", id="NP_000537.3", rettype="fasta")`, `efetch(db="snp",
id=["rs55863639", "rs587780067"])`). The complete list of parameters can be
found at ["The E-utilities In-Depth: Parameters, Syntax and
More"](https://www.ncbi.nlm.nih.gov/books/NBK25499/).

When a request succeeds the response object has a `data` field containing
formatted data, which can be saved to a file as demonstrated above. However,
users are often interested in a part of the response data and may want to
extract some fields in it. In such a case,
[EzXML.jl](https://github.com/bicycle1885/EzXML.jl) is helpful because it offers
lots of tools to handle XML documents. The first thing you need to do is
converting the response data into an XML document by `parsexml`:
```jlcon
julia> res = efetch(db="nuccore", id="NM_001126.3", retmode="xml")
Response(200 OK, 19 headers, 41536 bytes in body)

julia> using EzXML

julia> doc = parsexml(res.data)
EzXML.Document(EzXML.Node(<DOCUMENT_NODE@0x00007fdd4cc43770>))

```

After that, you can query fields you want using [XPath](https://en.wikipedia.org/wiki/XPath):
```jlcon
julia> seq = findfirst(doc, "/GBSet/GBSeq")
EzXML.Node(<ELEMENT_NODE@0x00007fdd49f34b10>)

julia> content(findfirst(seq, "GBSeq_definition"))
"Homo sapiens adenylosuccinate synthase (ADSS), mRNA"

julia> content(findfirst(seq, "GBSeq_accession-version"))
"NM_001126.3"

julia> length(find(seq, "//GBReference"))
10

julia> using Bio.Seq

julia> DNASequence(content(findfirst(seq, "GBSeq_sequence")))
2791nt DNA Sequence:
ACGGGAGTGGCGCGCCAGGCCGCGGAAGGGGCGTGGCCT…TGATTAAAAGAACCAAATATTTCTAGTATGAAAAAAAAA

```

Every function can take a context dictionary as its first argument to set
parameters into a query. Key-value pairs in a context are appended to a query in
addition to other parameters passed by keyword arguments.  The default context
is an empty dictionary that sets no parameters. This context dictionary is
especially useful when temporarily caching query UIDs into the Entrez History
server. A request to the Entrez system can be associated with cached data using
"WebEnv" and "query_key" parameters. In the following example, the search
results of `esearch` is saved in the Entrez History server (note
`usehistory=true`, which makes the server cache its search results) and then
their summaries are retrieved in the next call of `esummary`:
```jlcon
julia> context = Dict()  # create an empty context
Dict{Any,Any} with 0 entries

julia> res = esearch(context, db="pubmed", term="asthma[mesh] AND leukotrienes[mesh] AND 2009[pdat]", usehistory=true)
Response(200 OK, 18 headers, 1574 bytes in body)

julia> context  # the context dictionary has been updated
Dict{Any,Any} with 2 entries:
  :query_key => "1"
  :WebEnv    => "NCID_1_9251987_130.14.22.215_9001_1477389145_1960133…

julia> res = esummary(context, db="pubmed")  # retrieve summaries in context
Response(200 OK, 18 headers, 135463 bytes in body)

julia> write("asthma_leukotrienes_2009.xml", res.data)  # save data into a file
135463

shell> head asthma_leukotrienes_2009.xml
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1 20041029//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">
<eSummaryResult>
<DocSum>
        <Id>20113659</Id>
        <Item Name="PubDate" Type="Date">2009 Nov</Item>
        <Item Name="EPubDate" Type="Date"></Item>
        <Item Name="Source" Type="String">Zhongguo Dang Dai Er Ke Za Zhi</Item>
        <Item Name="AuthorList" Type="List">
                <Item Name="Author" Type="String">He MJ</Item>

```
