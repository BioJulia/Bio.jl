export Accession, EntrezGene, GenBank, RefSeq, GOTerm, browse


immutable Accession{T}
    data::T
end

=={T}(x::Accession{T}, y::Accession{T}) = x.data == y.data
=={T}(x::Accession{T}, y::String) = x.data == parse(T, y)
=={T}(x::String, y::Accession{T}) = parse(T, x) == y.data
hash(x::Accession) = hash(x.data)

parse{T}(::Type{Accession{T}}, x::String) = Accession(parse(T, x))
show(io::IO, x::Accession) = show(io, x.data)
browse(x::Accession) = browse(x.data)


# Entrez Gene

bitstype 32 EntrezGene

convert(::Type{EntrezGene}, geneid::Uint32) = box(EntrezGene, unbox(Uint32, geneid))
convert(::Type{Uint32}, geneid::EntrezGene) = box(Uint32, unbox(EntrezGene, geneid))

parse(::Type{EntrezGene}, s::String) = convert(EntrezGene, parse(Uint32, s))

function show(io::IO, geneid::EntrezGene)
    @printf io "%d" convert(Uint32, geneid)
end

hash(geneid::EntrezGene) = hash(convert(Uint32, geneid))

@osx_only function browse(geneid::EntrezGene)
    run(`open http://www.ncbi.nlm.nih.gov/gene/?term=$geneid`)
end


# GenBank

immutable GenBank
    accession::ASCIIString
    version::Uint8
end

==(x::GenBank, y::GenBank) = x.accession == y.accession && x.version == y.version

function parse(::Type{GenBank}, s::String)
    s = strip(s)
    dot = search(s, '.')
    @assert dot > 0
    acc = s[1:dot-1]
    version = parse(Uint8, s[dot+1:end])
    return GenBank(acc, version)
end

function show(io::IO, genbank::GenBank)
    print(io, genbank.accession)
    print(io, '.')
    print(io, genbank.version)
end

hash(genbank::GenBank) = hash(genbank.accession) $ hash(genbank.version)

@osx_only function browse(genbank::GenBank)
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

function parse(::Type{RefSeq}, s::String)
    s = strip(s)
    prefix = refseq_prefix_encode[s[1:3]]
    dot = search(s, '.')
    @assert dot > 0
    number = parse(Uint32, s[4:dot-1])
    version = parse(Uint8, s[dot+1:end])
    return RefSeq(prefix, number, version)
end

function show(io::IO, refseq::RefSeq)
    prefix = refseq_prefix_decode[refseq.prefix]
    nd = ndigits(refseq.number)
    if nd <= 6
        @printf io "%s%06d.%d" prefix refseq.number refseq.version
    elseif nd <= 9
        @printf io "%s%09d.%d" prefix refseq.number refseq.version
    else
        error("overflow")
    end
end

hash(refseq::RefSeq) = hash(refseq.prefix) $ hash(refseq.number) $ hash(refseq.version)

@osx_only function browse(refseq::RefSeq)
    run(`open http://www.ncbi.nlm.nih.gov/nuccore/$refseq`)
end


# Gene Ontology

bitstype 32 GOTerm

convert(::Type{GOTerm}, goterm::Uint32) = box(GOTerm, unbox(Uint32, goterm))
convert(::Type{Uint32}, goterm::GOTerm) = box(Uint32, unbox(GOTerm, goterm))

function parse(::Type{GOTerm}, s::String)
    s = strip(s)
    @assert startswith(s, "GO:")
    convert(GOTerm, parse(Uint32, s[4:end]))
end

function show(io::IO, goterm::GOTerm)
    @printf io "GO:%07d" convert(Uint32, goterm)
end

hash(goterm::GOTerm) = hash(convert(Uint32, goterm))

@osx_only function browse(goterm::GOTerm)
    run(`open http://amigo.geneontology.org/amigo/term/$goterm`)
end
