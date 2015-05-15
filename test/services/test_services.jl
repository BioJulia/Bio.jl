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
        @fact parse(Accession{sym}, "  $s \t") => accession
        @fact parse(Accession{sym}, convert(ASCIIString, s)) => accession
        @fact parse(Accession{sym}, convert(UTF8String, s)) => accession
        if guess
            @fact Accession(s) => accession
        end
    end

    function test_order(sym::Symbol, ss::Vector)
        accessions = [parse(Accession{sym}, s) for s in ss]
        for (s, accession) in zip(ss, sort(accessions))
            @fact string(accession) => s
        end
    end

    context("Entrez Gene") do
        list = split(chomp("""
        1
        912
        6598
        23681054
        105352660
        """))
        for s in list
            test_base(:EntrezGene, s)
        end
        test_order(:EntrezGene, list)
        @fact_throws parse(Accession{:EntrezGene}, "-1234") "negative"
        @fact_throws parse(Accession{:EntrezGene}, "01234") "starting with 0"
        @fact_throws parse(Accession{:EntrezGene}, "0x1234") "alphabet"
        @fact_throws parse(Accession{:EntrezGene}, "4294967296") "too large"
    end

    context("GenBank") do
        list = split(chomp("""
        DL128137
        DL128137.1
        JL971710
        JL971710.1
        JL971710.2
        U46667.1
        """))
        for s in list
            test_base(:GenBank, s)
        end
        test_order(:GenBank, list)
        @fact_throws parse(Accession{:GenBank}, "12345") "no prefix"
    end

    context("RefSeq") do
        list = split(chomp("""
        NC_000022.11
        NG_042145.1
        NW_012162423.1
        NZ_LCTL00000000
        NZ_LCTL01001056.1
        XM_011530345.1
        XP_011528647
        XP_011528647.1
        """))
        for s in list
            test_base(:RefSeq, s, guess=true)
        end
        test_order(:RefSeq, list)
        @fact_throws parse(Accession{:RefSeq}, "XX_000022.1") "invalid prefix"
    end

    context("Consensus CDS (CCDS)") do
        list = split(chomp("""
        CCDS1.1
        CCDS5251
        CCDS5251.1
        CCDS5251.2
        """))
        for s  in list
            test_base(:CCDS, s, guess=true)
        end
        test_order(:CCDS, list)
        @fact_throws parse(Accession{:CCDS}, "CDS1234.1") "invalid prefix"
    end

    context("Ensembl") do
        list = split(chomp("""
        ENSG00000139618
        ENSG00000139618.13
        """))
        for s in list
            test_base(:Ensembl, s, guess=true)
        end
        test_order(:Ensembl, list)
    end

    context("Gene Ontology") do
        list = split(chomp("""
        GO:0000003
        GO:0016514
        GO:0044183
        GO:0044848
        GO:2001306
        """))
        for s in list
            test_base(:GeneOntology, s, guess=true)
        end
        test_order(:GeneOntology, list)
        @fact_throws parse(Accession{:GeneOntology}, "GO:123456") "short number"
        @fact_throws parse(Accession{:GeneOntology}, "GO:12345678") "long number"
        @fact_throws parse(Accession{:GeneOntology}, "GX:1234567") "invalid prefix"
    end

    context("UniProt") do
        list = split(chomp("""
        A0A022YWF9
        A2BC19
        N1QTE2
        O55743
        P12345
        Q197B5
        Q6GZX4
        """))
        for s in list
            test_base(:UniProt, s, guess=true)
        end
        test_order(:UniProt, list)
        @fact_throws parse(Accession{:UniProt}, "P1234") "short number"
        @fact_throws parse(Accession{:UniProt}, "P123456") "long number"
    end
end

end # TestServices
