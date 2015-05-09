module TestServices

using FactCheck
using Bio
using Bio.Services

facts("Accession") do
    context("Entrez Gene") do
        # SMARCB1: Homo sapiens
        smarcb1 = parse(Accession{EntrezGene}, "6598")
        @fact string(smarcb1) => "6598"
        @fact smarcb1 == "6598" => true
        @fact "6598" == smarcb1 => true
        # SMARCB1: Zonotrichia albicollis
        smarcb1 = parse(Accession{EntrezGene}, "102067278")
        @fact string(smarcb1) => "102067278"
    end

    context("GenBank") do
        x = parse(Accession{GenBank}, "U46667.1")
        @fact x == "U46667.1" => true
        @fact string(x) => "U46667.1"
    end

    context("RefSeq") do
        @fact parse(Accession{RefSeq}, "NC_000022.11") |> string => "NC_000022.11"
        @fact parse(Accession{RefSeq}, "XM_011530345.1") |> string => "XM_011530345.1"
        @fact parse(Accession{RefSeq}, "XP_011528647.1") |> string => "XP_011528647.1"
    end

    context("Gene Ontology") do
        # SWI/SNF complex
        @fact parse(Accession{GOTerm}, "GO:0016514") |> string => "GO:0016514"
    end
end

end # TestServices
