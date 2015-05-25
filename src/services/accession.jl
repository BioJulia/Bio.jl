# Accessions
# ==========

export Accession, browse, uri

using URIParser

# Accession number type.
#   S: Symbol    database name
#   T: DataType  encoding type of accession number
immutable Accession{S,T}
    # encoded accession number
    data::T
end

Accession(s::String) = parse(Accession, s)

if VERSION >= v"v0.4-"
    Base.call{S}(::Type{Accession{S}}, s::String) = parse(Accession{S}, s)
end

=={S,T}(x::Accession{S,T}, y::Accession{S,T}) = x.data == y.data
=={S}(x::Accession{S}, y::String) = x == parse(Accession{S}, y)
=={S}(x::String, y::Accession{S}) = parse(Accession{S}, x) == y
isless{S,T}(x::Accession{S,T}, y::Accession{S,T}) = isless(x.data, y.data)
hash(x::Accession) = hash(x.data)

function browse(accession::Accession)
    link = uri(accession, format=:browser)
    # copied from Gadfly.jl
    if OS_NAME == :Darwin
        run(`open $link`)
    elseif OS_NAME == :Linux || OS_NAME == :FreeBSD
        run(`xdg-open $link`)
    elseif OS_NAME == :Windows
        run(`$(ENV["COMSPEC"]) /c start $link`)
    else
        warn("Opening a web browser is not supported on OS $(string(OS_NAME))")
    end
end

# Smart Accession Parser
#
# This function should not be used outside of an interactive session because
# the strategy of guessing accession types is ad hoc and the results should be
# checked by callers: the inferred types may not be what you thought!
function parse(::Type{Accession}, s::String)
    if ismatch(Accession{:RefSeq}, s)
        return parse(Accession{:RefSeq}, s)
    elseif ismatch(Accession{:CCDS}, s)
        return parse(Accession{:CCDS}, s)
    elseif ismatch(Accession{:Ensembl}, s)
        return parse(Accession{:Ensembl}, s)
    elseif ismatch(Accession{:GeneOntology}, s)
        return parse(Accession{:GeneOntology}, s)
    elseif ismatch(Accession{:UniProt}, s)
        return parse(Accession{:UniProt}, s)
    end
    error("cannot guess accession number type")
end

macro check_match(name, s)
    quote
        if !ismatch(Accession{$name}, $s)
            error("invalid $($name) accession number: '$($s)'")
        end
    end
end


# Accession Encodings
# ===================

# Accession numbers must support following methods:
#   * ==(x::T, y::T)
#   * isless(x::T, y::T)
#   * hash(x::T)
#   * ismatch(::Type{Accession{S}}, s::String)
#   * parse(::Type{Accession{S}}, s::String)
#   * show(io::IO, accnum::Accession{S})
#
# And recommended methods:
#   * uri(accnum::Accession{S}; format::Symbol=:browser)


# internal representation of versioned accession numbers
immutable Versioned{T}
    accession::T
    version::Nullable{Uint8}
end

Versioned(x) = Versioned(x, Nullable{Uint8}())
Versioned(x, v) = Versioned(x, parse(Uint8, v))

function =={T}(x::Versioned{T}, y::Versioned{T})
    if x.accession != y.accession
        return false
    elseif isnull(x.version)
        return ifelse(isnull(y.version), true, false)
    elseif isnull(y.version)
        return false
    end
    return get(x.version) == get(y.version)
end

function isless{T}(x::Versioned{T}, y::Versioned{T})
    if x.accession != y.accession
        return x.accession < y.accession
    elseif isnull(x.version)
        return ifelse(isnull(y.version), false, true)
    elseif isnull(y.version)
        return false
    end
    return get(x.version) < get(y.version)
end

hash(x::Versioned) = hash(x.accession) $ hash(x.version)

isversioned(x::Versioned) = !isnull(x.version)
getversion(x::Versioned) = get(x.version)


# Entrez Gene
# -----------

# webpage:
#   http://www.ncbi.nlm.nih.gov/gene
# accession format:
#   http://www.ncbi.nlm.nih.gov/books/NBK3841/#EntrezGene.Numbering_system

function ismatch(::Type{Accession{:EntrezGene}}, s::String)
    # typemax(UInt32) == 4,294,967,295
    return ismatch(r"^\s*(:?[1-9]\d{0,8}|[1-3]?\d{9})\s*$", s)
end

function parse(::Type{Accession{:EntrezGene}}, s::String)
    @check_match :EntrezGene s
    n = parse(Uint32, s)
    return Accession{:EntrezGene,Uint32}(n)
end

function show(io::IO, geneid::Accession{:EntrezGene})
    @printf io "%d" convert(Uint32, geneid.data)
end

function uri(geneid::Accession{:EntrezGene}; format::Symbol=:browser)
    @assert format === :browser
    return URI("http://www.ncbi.nlm.nih.gov/gene/?term=$geneid")
end


# GenInfo Identifier (GI)
# -----------------------

# accession format:
#   http://www.ncbi.nlm.nih.gov/genbank/sequenceids/

function ismatch(::Type{Accession{:GI}}, s::String)
    return ismatch(r"^\s*[1-9]\d*\s*$", s)
end

function parse(::Type{Accession{:GI}}, s::String)
    @check_match :GI s
    return Accession{:GI,Uint}(parse(Uint, s))
end

show(io::IO, gi::Accession{:GI}) = print(io, gi.data)

function uri(gi::Accession{:GI}; format::Symbol=:browser)
    @assert format === :browser
    return URI("http://www.ncbi.nlm.nih.gov/nuccore/$gi")
end


# GenBank 
# -------

# webpage:
#   http://www.ncbi.nlm.nih.gov/genbank/
# accession format:
#   http://www.ncbi.nlm.nih.gov/Sequin/acc.html

function ismatch(::Type{Accession{:GenBank}}, s::String)
    return ismatch(r"^\s*[A-Z]{1,5}\d{5,10}(:?\.\d+)?\s*$", s)
    #                    ^~~~~~~~~~^~~~~~~~^~~~~~~~~
    #                    prefix    |       version
    #                              numerals
end

function parse(::Type{Accession{:GenBank}}, s::String)
    @check_match :GenBank s
    dot = search(s, '.')
    if dot == 0
        accession = strip(s)
        version = Nullable{Uint8}()
    else
        accession = strip(s[1:dot-1])
        version = Nullable(parse(Uint8, s[dot+1:end]))
    end
    return Accession{:GenBank,Versioned{ASCIIString}}(Versioned{ASCIIString}(accession, version))
end

function show(io::IO, genbank::Accession{:GenBank})
    write(io, genbank.data.accession)
    if isversioned(genbank.data)
        write(io, '.')
        write(io, string(getversion(genbank.data)))
    end
    return
end

function uri(genbank::Accession{:GenBank}; format::Symbol=:browser)
    @assert format === :browser
    return URI("http://www.ncbi.nlm.nih.gov/nuccore/$genbank")
end


# RefSeq
# ------

# webpage:
#   http://www.ncbi.nlm.nih.gov/refseq/
# accession format:
#   http://www.ncbi.nlm.nih.gov/books/NBK21091/
#   http://www.ncbi.nlm.nih.gov/genbank/sequenceids/

function ismatch(::Type{Accession{:RefSeq}}, s::String)
    return ismatch(r"^\s*(:?A[CP]|N[CGTWSZMRP]|X[MRP]|YP|ZP)_[A-Z0-9]+(:?\.\d+)?\s*$", s)
    #                    ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^~~~~~~~~^~~~~~~~~~
    #                    prefix                              |        version (optional)
    #                                                        alpha-numerals
end

function parse(::Type{Accession{:RefSeq}}, s::String)
    @check_match :RefSeq s
    dot = search(s, '.')
    if dot == 0
        accession = strip(s)
        version = Nullable{Uint8}()
    else
        accession = strip(s[1:dot-1])
        version = Nullable(parse(Uint8, s[dot+1:end]))
    end
    return Accession{:RefSeq,Versioned{ASCIIString}}(Versioned{ASCIIString}(accession, version))
end

function show(io::IO, refseq::Accession{:RefSeq})
    write(io, refseq.data.accession)
    if isversioned(refseq.data)
        write(io, '.')
        write(io, string(getversion(refseq.data)))
    end
    return
end

function uri(refseq::Accession{:RefSeq}; format::Symbol=:browser)
    @assert format === :browser
    return URI("http://www.ncbi.nlm.nih.gov/nuccore/$refseq")
end


# Consensus CDS (CCDS)
# --------------------

# webpage:
#   http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi
# accession format:
#   http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi#ccdsIds

function ismatch(::Type{Accession{:CCDS}}, s::String)
    return ismatch(r"^\s*CCDS\d+(:?\.\d+)?\s*$", s)
end

function parse(::Type{Accession{:CCDS}}, s::String)
    @check_match :CCDS s
    ccds = search(s, "CCDS")
    dot = search(s, '.')
    if dot == 0
        accession = parse(Uint32, s[last(ccds)+1:end])
        version = Nullable{Uint8}()
    else
        accession = parse(Uint32, s[last(ccds)+1:dot-1])
        version = Nullable(parse(Uint8, s[dot+1:end]))
    end
    return Accession{:CCDS,Versioned{Uint32}}(Versioned{Uint32}(accession, version))
end

function show(io::IO, ccds::Accession{:CCDS})
    @printf io "CCDS%d" ccds.data.accession
    if isversioned(ccds.data)
        write(io, '.')
        write(io, string(getversion(ccds.data)))
    end
end

function uri(ccds::Accession{:CCDS}; format::Symbol=:browser)
    @assert format === :browser
    return URI("http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=$ccds")
end


# Ensembl
# -------

# webpage:
#   http://www.ensembl.org/index.html
# accession format:
#   http://www.ensembl.org/info/genome/stable_ids/index.html

function ismatch(::Type{Accession{:Ensembl}}, s::String)
    return ismatch(r"^\s*(:?ENS(:?[A-Z]{3,4}?)?|FB)(:?E|FM|G|GT|P|R|T)\d+(:?\.\d+)?\s*$", s)
    #                    ^~~~~~~~~~~~~~~~~~~~~~~~~~^~~~~~~~~~~~~~~~~~~^~~^~~~~~~~~
    #                    species prefix            feature prefix     |  version (optional)
    #                                                                 numerals
end

function parse(::Type{Accession{:Ensembl}}, s::String)
    @check_match :Ensembl s
    dot = search(s, '.')
    if dot == 0
        accession = strip(s)
        version = Nullable{Uint8}()
    else
        accession = strip(s[1:dot-1])
        version = Nullable(parse(Uint8, s[dot+1:end]))
    end
    return Accession{:Ensembl,Versioned{ASCIIString}}(Versioned{ASCIIString}(accession, version))
end

function show(io::IO, ensembl::Accession{:Ensembl})
    write(io, ensembl.data.accession)
    if isversioned(ensembl.data)
        write(io, '.')
        write(io, string(getversion(ensembl.data)))
    end
    return
end

# TODO: Is there a permanent link to an entry?
#function uri(ensembl::Accession{:Ensembl}; format::Symbol=:browser)
#end


# Gene Ontology
# -------------

# webpage:
#   http://geneontology.org/
# accession format:
#   http://geneontology.org/page/ontology-structure#termstru

function ismatch(::Type{Accession{:GeneOntology}}, s::String)
    return ismatch(r"^\s*GO:\d{7}\s*$", s)
end

function parse(::Type{Accession{:GeneOntology}}, s::String)
    @check_match :GeneOntology s
    go = search(s, "GO:")
    n = parse(Uint32, s[last(go)+1:end])
    return Accession{:GeneOntology,Uint32}(n)
end

function show(io::IO, go::Accession{:GeneOntology})
    @printf io "GO:%07d" convert(Uint32, go.data)
end

function uri(go::Accession{:GeneOntology}; format::Symbol=:browser)
    if format === :browser
        q = "id=$go"
    else
        @assert format ∈ [:mini, :obo, :oboxml]
        q = "id=$go&format=$format"
    end
    # TODO: AmiGO or QuickGO
    return URI("http://www.ebi.ac.uk/QuickGO/GTerm?$q")
end


# UniProt
# -------

# webpage:
#   http://www.uniprot.org/
# accession format:
#   http://www.uniprot.org/help/accession_numbers

function ismatch(::Type{Accession{:UniProt}}, s::String)
    return ismatch(r"^\s*(:?[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](:?[A-Z][A-Z0-9]{2}[0-9]){1,2})\s*$", s)
end

function parse(::Type{Accession{:UniProt}}, s::String)
    @check_match :UniProt s
    return Accession{:UniProt,ASCIIString}(strip(s))
end

function show(io::IO, uniprot::Accession{:UniProt})
    write(io, uniprot.data)
    return
end

function uri(uniprot::Accession{:UniProt}; format::Symbol=:browser)
    if format === :browser
        ext = ""
    else
        @assert format ∈ [:txt, :fasta, :xml, :rdf, :gff]
        ext = ".$format"
    end
    return URI("http://www.uniprot.org/uniprot/$uniprot$ext")
end
