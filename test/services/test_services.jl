module TestServices

using FactCheck
using Bio
using Bio.Services

facts("Accession Number") do
    function test_base(sym::Symbol, s::String; guess=false)
        accession = parse(Accession{sym}, s)
        @fact isa(accession, Accession{sym}) => true
        @fact string(accession) => strip(s)
        @fact accession == s => true
        @fact s == accession => true
        @fact parse(Accession{sym}, s) => accession
        @fact parse(Accession{sym}, convert(ASCIIString, s)) => accession
        @fact parse(Accession{sym}, convert(UTF8String, s)) => accession
        if guess
            @fact Accession(s) => accession
        end
    end

    context("Entrez Gene") do
        list = """
        1
        6598
        105352660
         6598 \t 
        """ |> chomp
        for s in split(list, '\n')
            test_base(:EntrezGene, s)
        end
        @fact_throws parse(Accession{:EntrezGene}, "-1234") "negative"
        @fact_throws parse(Accession{:EntrezGene}, "0x1234") "alphabet"
        @fact_throws parse(Accession{:EntrezGene}, "4294967296") "too large"
    end

    context("GenBank") do
        list = """
        U46667.1
        DL128137.1
         DL128137.1 \t 
        DL128137
        """ |> chomp
        for s in split(list, '\n')
            test_base(:GenBank, s)
        end
    end

    context("RefSeq") do
        list = """
        NC_000022.11
        XM_011530345.1
        XP_011528647.1
         XP_011528647.1 \t 
        XP_011528647
        """ |> chomp
        for s in split(list, '\n')
            test_base(:RefSeq, s, guess=true)
        end
        @fact_throws parse(Accession{:RefSeq}, "XX_000022.1") "invalid prefix"
    end

    context("Consensus CDS (CCDS)") do
        list = """
        CCDS1.1
        CCDS5251.2
         CCDS5251.2 \t 
        CCDS5251
        """ |> chomp
        for s  in split(list, '\n')
            test_base(:CCDS, s, guess=true)
        end
        @fact_throws parse(Accession{:CCDS}, "CDS1234.1") "invalid prefix"
    end

    context("Ensembl") do
        list = """
        ENSG00000139618.13
         ENSG00000139618.13 \t 
        ENSG00000139618
        """ |> chomp
        for s in split(list, '\n')
            test_base(:Ensembl, s, guess=true)
        end
    end

    context("Gene Ontology") do
        list = """
        GO:0016514
        GO:0044848
         GO:0016514 \t 
        """ |> chomp
        for s in split(list, '\n')
            test_base(:GeneOntology, s, guess=true)
        end
        @fact_throws parse(Accession{:GeneOntology}, "GO:123456") "short number"
        @fact_throws parse(Accession{:GeneOntology}, "GO:12345678") "long number"
        @fact_throws parse(Accession{:GeneOntology}, "GX:1234567") "invalid prefix"
    end

    context("UniProt") do
        list = """
        A2BC19
        P12345
        A0A022YWF9
         A0A022YWF9 \t 
        """ |> chomp
        for s in split(list, '\n')
            test_base(:UniProt, s, guess=true)
        end
        @fact_throws parse(Accession{:UniProt}, "P1234") "short number"
        @fact_throws parse(Accession{:UniProt}, "P123456") "long number"
    end
end

end # TestServices
