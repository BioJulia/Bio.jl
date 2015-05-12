export Accession, browse

# Accession numbers (S::Symbol database name, T::DataType encoding type)
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
hash(x::Accession) = hash(x.data)

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
    end
    error("cannot guess accession number type")
end

# parser of unsigne decimal integer; this is very common among accession numbers
function parse_decimal_uint(s::String, start::Int=1, stop::Int=endof(s))
    @assert 1 <= start <= stop <= endof(s)
    n = zero(Uint)
    for i in start:stop
        c = s[i]
        if '0' <= c <= '9'
            n′ = 10n + (c - '0')
            if n′ < n
                throw(OverflowError())
            end
            n = n′
        elseif isspace(c)
            break
        else
            error("invalid character: '$c'")
        end
    end
    return n
end

function search_nonspace(s::String)
    for i in 1:endof(s)
        if !isspace(s[i])
            return i
        end
    end
    return 0
end


# Entrez Gene

function parse(::Type{Accession{:EntrezGene}}, s::String)
    n = parse_decimal_uint(s, search_nonspace(s))
    if n > typemax(Uint32)
        error("too large accession number")
    end
    return Accession{:EntrezGene,Uint32}(n)
end

function show(io::IO, geneid::Accession{:EntrezGene})
    @printf io "%d" convert(Uint32, geneid.data)
end

@osx_only function browse(geneid::Accession{:EntrezGene})
    run(`open http://www.ncbi.nlm.nih.gov/gene/?term=$geneid`)
end


# GenBank

# http://www.ncbi.nlm.nih.gov/Sequin/acc.html

immutable GenBank
    accession::ASCIIString
    version::Uint8
end

==(x::GenBank, y::GenBank) = x.version == y.version && x.accession == y.accession

function parse(::Type{Accession{:GenBank}}, s::String)
    s = strip(s)
    dot = search(s, '.')
    @assert dot > 0
    accession = s[1:dot-1]
    version = parse_decimal_uint(SubString(s, dot+1))
    return Accession{:GenBank,GenBank}(GenBank(accession, version))
end

function show(io::IO, genbank::Accession{:GenBank})
    print(io, genbank.data.accession)
    print(io, '.')
    print(io, genbank.data.version)
end

hash(genbank::GenBank) = hash(genbank.accession) $ hash(genbank.version)

@osx_only function browse(genbank::Accession{:GenBank})
    run(`open http://www.ncbi.nlm.nih.gov/nuccore/$genbank`)
end


# RefSeq

immutable RefSeq
    prefix::Uint8
    number::Uint32
    version::Uint8
end

# generate lookup table
macro gen_table(name, typ, xs)
    xs = eval(xs)
    enc_table = symbol(string(name, :_encode))
    dec_table = symbol(string(name, :_decode))
    esc(quote
        const $enc_table = [x => convert($typ, i) for (i, x) in enumerate($xs)]
        const $dec_table = $xs
    end)
end

# http://www.ncbi.nlm.nih.gov/books/NBK21091/ (Table 1)
@gen_table refseq_prefix Uint8 [
    "AC_", "NC_", "NG_", "NT_", "NW_",
    "NS_", "NZ_", "NM_", "NR_", "XM_",
    "XR_", "AP_", "NP_", "YP_", "XP_",
    "ZP_",
]

function ismatch(::Type{Accession{:RefSeq}}, s::String)
    return ismatch(r"^\s*(:?A[CP]|N[CGTWSZMRP]|X[MRP]|YP|ZP)_\d{6}(:?\d{3})?\.\d+\s*$", s)
    #                    ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^~~~~~~~~~~~~~~   ^~
    #                    prefix                              number            version
end

function parse(::Type{Accession{:RefSeq}}, s::String)
    s = strip(s)
    prefix = refseq_prefix_encode[s[1:3]]
    dot = search(s, '.')
    @assert dot > 0
    number = parse(Uint32, s[4:dot-1])
    version = parse(Uint8, s[dot+1:end])
    return Accession{:RefSeq,RefSeq}(RefSeq(prefix, number, version))
end

function show(io::IO, refseq::Accession{:RefSeq})
    prefix = refseq_prefix_decode[refseq.data.prefix]
    nd = ndigits(refseq.data.number)
    if nd <= 6
        @printf io "%s%06d.%d" prefix refseq.data.number refseq.data.version
    elseif nd <= 9
        @printf io "%s%09d.%d" prefix refseq.data.number refseq.data.version
    else
        error("number overflow")
    end
end

hash(refseq::RefSeq) = hash(refseq.prefix) $ hash(refseq.number) $ hash(refseq.version)

@osx_only function browse(refseq::Accession{:RefSeq})
    run(`open http://www.ncbi.nlm.nih.gov/nuccore/$refseq`)
end

# Consensus CDS (CCDS)

# http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi#ccdsIds

immutable CCDS
    number::Uint32
    version::Uint8
end

function ismatch(::Type{Accession{:CCDS}}, s::String)
    return ismatch(r"^\s*CCDS\d+\.\d+\s*$", s)
end


# Ensembl

# http://www.ensembl.org/info/genome/stable_ids/index.html

function ismatch(::Type{Accession{:Ensembl}}, s::String)
    return ismatch(r"^\s*(:?ENS(:?[A-Z]{3,4}?)|FB)(:?E|FM|G|GT|P|R|T)\d+\s*$", s)
    #                    ^~~~~~~~~~~~~~~~~~~~~~~~~^~~~~~~~~~~~~~~~~~~^~~
    #                    species prefix           feature prefix     number
end

function parse(::Type{Accession{:Ensembl,ASCIIString}}, s::String)
    return Accession{:Ensembl,ASCIIString}(strip(s))
end


# Gene Ontology

function ismatch(::Type{Accession{:GeneOntology}}, s::String)
    return ismatch(r"^\s*GO:\d{7}\s*$", s)
end

function parse(::Type{Accession{:GeneOntology}}, s::String)
    if !ismatch(Accession{:GeneOntology}, s)
        error("invalid GeneOntology accession number")
    end
    i = search_nonspace(s)
    # `+ 3` is the offset of the prefix 'GO:'
    n = parse_decimal_uint(s, i + 3)
    return Accession{:GeneOntology,Uint32}(n)
end

function show(io::IO, go::Accession{:GeneOntology})
    @printf io "GO:%07d" convert(Uint32, go.data)
end

@osx_only function browse(go::Accession{:GeneOntology})
    run(`open http://amigo.geneontology.org/amigo/term/$go`)
end
