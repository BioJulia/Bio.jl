module TestServices

using FactCheck
using Bio
using Bio.Services

facts("Accession Number") do
    context("Entrez Gene") do
        list = """
        1
        6598
        105352660
         6598 \t 
        """ |> chomp
        for s in split(list, '\n')
            geneid = parse(Accession{:EntrezGene}, s)
            @fact isa(geneid, Accession{:EntrezGene}) => true
            @fact string(geneid) => strip(s)
            @fact geneid == s => true
            @fact s == geneid => true
            @fact geneid == parse(Accession{:EntrezGene}, s) => true
        end
        @fact_throws parse(Accession{:EntrezGene}, "-1234") "negative"
        @fact_throws parse(Accession{:EntrezGene}, "0x1234") "alphabet"
        # TODO: this fails
        @fact_throws parse(Accession{:EntrezGene}, "4294967296") "too large"
    end

    context("GenBank") do
        list = """
        U46667.1
        DL128137.1
         DL128137.1 \t 
        """ |> chomp
        for s in split(list, '\n')
            genbank = parse(Accession{:GenBank}, s)
            @fact isa(genbank, Accession{:GenBank}) => true
            @fact string(genbank) => strip(s)
            @fact genbank == s => true
            @fact s == genbank => true
            @fact genbank == parse(Accession{:GenBank}, s) => true
        end
    end

    context("RefSeq") do
        list = """
        NC_000022.11
        XM_011530345.1
        XP_011528647.1
         XP_011528647.1 \t 
        """ |> chomp
        for s in split(list, '\n')
            refseq = parse(Accession{:RefSeq}, s)
            @fact isa(refseq, Accession{:RefSeq}) => true
            @fact string(refseq) => strip(s)
            @fact refseq == s => true
            @fact s == refseq => true
            @fact refseq == parse(Accession{:RefSeq}, s) => true
            @fact refseq == Accession(s) => true
        end
    end

    context("Gene Ontology") do
        list = """
        GO:0016514
        GO:0044848
         GO:0016514 \t 
        """ |> chomp
        for s in split(list, '\n')
            go = parse(Accession{:GeneOntology}, s)
            @fact isa(go, Accession{:GeneOntology}) => true
            @fact string(go) => strip(s)
            @fact go == s => true
            @fact s == go => true
            @fact go == parse(Accession{:GeneOntology}, s) => true
            @fact go == Accession(s) => true
        end
    end
end

end # TestServices
