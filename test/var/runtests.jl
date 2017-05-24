module TestVar

using Base.Test

using Bio: Seq, Var
using TestFunctions
using PairwiseListMatrices
using IntervalTrees: IntervalValue
import BufferedStreams: BufferedInputStream
import YAML

function generate_possibilities_tester{A<:Union{DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4}, RNAAlphabet{2}}}(::Type{A})
    symbols = alphabet(A)
    arra = Vector{eltype(A)}()
    arrb = Vector{eltype(A)}()
    for i in 1:length(symbols), j in i:length(symbols)
        push!(arra, symbols[i])
        push!(arrb, symbols[j])
    end
    return BioSequence{A}(arra), BioSequence{A}(arrb)
end

for alph in (:DNAAlphabet, :RNAAlphabet)
    @eval function generate_possibilities_tester(::Type{$alph{2}}, ::Type{$alph{4}})
        arra = Vector{eltype($alph)}()
        arrb = Vector{eltype($alph)}()
        for i in alphabet($alph{2}), j in alphabet($alph{4})
            push!(arra, i)
            push!(arrb, j)
        end
        return BioSequence{$alph{2}}(arra), BioSequence{$alph{4}}(arrb)
    end
end

#=
@testset "Site counting" begin
    @testset "Naive methods" begin

        alphabets = (DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4}, RNAAlphabet{2})
        NC = Seq.NaiveCount

        for alph in alphabets

            # Answers to these tests were worked out manually to verify count_sites_naive was working.
            # seqA and seqB contain all possible observations of sites.

            istwobit = Seq.bitsof(alph) == 2

            seqA, seqB = generate_possibilities_tester(alph)

            # Test when sequences are of the same bitencoding.

            @test count(Mutated, NC, seqA, seqB) == count(Mutated, NC, seqB, seqA) == ifelse(istwobit, 6, (6, 10))
            @test count(Conserved, NC, seqA, seqB) == count(Conserved, NC, seqB, seqA) == ifelse(istwobit, 4, (4, 10))
            @test count(Transition, NC, seqA, seqB) == count(Transition, NC, seqB, seqA) == ifelse(istwobit, 2, (2, 10))
            @test count(Transversion, NC, seqA, seqB) == count(Transversion, NC, seqB, seqA) == ifelse(istwobit, 4, (4, 10))
        end

        # Test for when sequences are of different bitencodings.
        for alphs in [(DNAAlphabet{2}, DNAAlphabet{4}),
                      (RNAAlphabet{2}, RNAAlphabet{4})]
            seqA, seqB = generate_possibilities_tester(alphs...)

            @test count(Mutated, NC, seqA, seqB) == count(Mutated, NC, seqB, seqA) == (12, 16)
            @test count(Conserved, NC, seqA, seqB) == count(Conserved, NC, seqB, seqA) == (4, 16)
            @test count(Transition, NC, seqA, seqB) == count(Transition, NC, seqB, seqA) == (4, 16)
            @test count(Transversion, NC, seqA, seqB) == count(Transversion, NC, seqB, seqA) == (8, 16)
        end
    end

    @testset "Bit parallel methods" begin
        @testset "4 bit encoding" begin
            NC = Seq.NaiveCount
            BC = Seq.BitparCount
            alphabets = (DNAAlphabet{4}, RNAAlphabet{4})
            for alph in alphabets
                for _ in 1:50
                    seqA = random_seq(alph, rand(10:100))
                    seqB = random_seq(alph, rand(10:100))
                    subA = seqA[1:rand(10:length(seqA))]
                    subB = seqB[1:rand(10:length(seqB))]
                    @test count(Mutated, BC, subA, subB) == count(Mutated, BC, subB, subA) == count(Mutated, NC, subA, subB)
                    @test count(Conserved, BC, subA, subB) == count(Conserved, BC, subB, subA) == count(Conserved, NC, subA, subB)
                end
            end
        end

        @testset "2 bit encoding" begin
            NC = Seq.NaiveCount
            BC = Seq.BitparCount
            alphabets = (DNAAlphabet{2}, RNAAlphabet{2})
            for alph in alphabets
                for _ in 1:50
                    seqA = random_seq(alph, rand(10:100))
                    seqB = random_seq(alph, rand(10:100))
                    subA = seqA[1:rand(10:length(seqA))]
                    subB = seqB[1:rand(10:length(seqB))]
                    @test count(Mutated, BC, subA, subB) == count(Mutated, BC, subB, subA) == count(Mutated, NC, subA, subB)
                    @test count(Conserved, BC, subA, subB) == count(Conserved, BC, subB, subA) == count(Conserved, NC, subA, subB)
                end
            end
        end

    end
        


        @testset "Pairwise methods" begin
            dnas = [dna"ATCGCCA-",
                    dna"ATCGCCTA",
                    dna"ATCGCCT-",
                    dna"GTCGCCTA"]

            rnas = [rna"AUCGCCA-",
                    rna"AUCGCCUA",
                    rna"AUCGCCU-",
                    rna"GUCGCCUA"]

            for i in (dnas, rnas)
                @test count_pairwise(CONSERVED, i...) == PairwiseListMatrix([(6,7), (6,7), (5,7), (7,7), (7,8), (6,7)], false)

                @test count_pairwise(MUTATED, i...) == PairwiseListMatrix([(1,7), (1,7), (2,7), (0,7), (1,8), (1,7)], false)

                @test count_pairwise(TRANSITION, i...) == PairwiseListMatrix([(0,7), (0,7), (1,7), (0,7), (1,8), (1,7)], false)

                @test count_pairwise(TRANSVERSION, i...) == PairwiseListMatrix([(1,7), (1,7), (1,7), (0,7), (0,8), (0,7)], false)
            end
        end

        @testset "Windowed methods" begin
            dnaA = dna"ATCGCCA-M"
            dnaB = dna"ATCGCCTAA"
            rnaA = rna"AUCGCCA-M"
            rnaB = rna"AUCGCCUAA"
            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count(CONSERVED, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (3,3)),
                                                                   IntervalValue(2, 4, (3,3)),
                                                                   IntervalValue(3, 5, (3,3)),
                                                                   IntervalValue(4, 6, (3,3)),
                                                                   IntervalValue(5, 7, (2,3)),
                                                                   IntervalValue(6, 8, (1,2)),
                                                                   IntervalValue(7, 9, (0,1))]
                @test count(MUTATED, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (0,3)),
                                                                 IntervalValue(2, 4, (0,3)),
                                                                 IntervalValue(3, 5, (0,3)),
                                                                 IntervalValue(4, 6, (0,3)),
                                                                 IntervalValue(5, 7, (1,3)),
                                                                 IntervalValue(6, 8, (1,2)),
                                                                 IntervalValue(7, 9, (1,1))]
                @test count(TRANSITION, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (0,3)),
                                                                    IntervalValue(2, 4, (0,3)),
                                                                    IntervalValue(3, 5, (0,3)),
                                                                    IntervalValue(4, 6, (0,3)),
                                                                    IntervalValue(5, 7, (0,3)),
                                                                    IntervalValue(6, 8, (0,2)),
                                                                    IntervalValue(7, 9, (0,1))]
                @test count(TRANSVERSION, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (0,3)),
                                                                      IntervalValue(2, 4, (0,3)),
                                                                      IntervalValue(3, 5, (0,3)),
                                                                      IntervalValue(4, 6, (0,3)),
                                                                      IntervalValue(5, 7, (1,3)),
                                                                      IntervalValue(6, 8, (1,2)),
                                                                      IntervalValue(7, 9, (1,1))]
            end
        end

        @testset "MASH distances" begin
            a = minhash(dna"ATCGCCA-", 4, 3)
            b = minhash(dna"ATCGCCTA", 4, 3)
            @test_approx_eq_eps mashdistance(a, b) 0.2745 1e-3
            @test mashdistance(a, a) == 0
            @test a.sketch == sort(a.sketch)
        end
    end

end
=#
#=
@testset "Distance Computation" begin

    dnas1 = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
    m1 = seqmatrix(dnas1, :seq)

    dnas2 = [dna"attgaacctggntttccgaa",
             dna"atacagagtatacrgtcgtc"]
    dnas3 = [dna"attgaacctgtntttccgaa",
             dna"atagaacgtatatrgccgtc"]
    m2 = seqmatrix(dnas2, :seq)

    @test distance(Count{AnyMutation}, dnas1) == ([12], [16])
    @test distance(Count{TransitionMutation}, dnas1) == ([4], [16])
    @test distance(Count{TransversionMutation}, dnas1) == ([8], [16])
    @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])
    @test distance(Count{AnyMutation}, m1) == ([12], [16])
    @test distance(Count{TransitionMutation}, m1) == ([4], [16])
    @test distance(Count{TransversionMutation}, m1) == ([8], [16])
    @test distance(Count{Kimura80}, m1) == ([4], [8], [16])

    @test distance(Count{AnyMutation}, dnas2, 5, 5)[1][:] == [2, 4, 3, 3]
    @test distance(Count{AnyMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
    @test distance(Count{TransitionMutation}, dnas2, 5, 5)[1][:] == [0, 2, 1, 1]
    @test distance(Count{TransitionMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
    @test distance(Count{TransversionMutation}, dnas2, 5, 5)[1][:] == [2, 2, 2, 2]
    @test distance(Count{TransversionMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
    @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])

    @test distance(Count{AnyMutation}, dnas2) == ([12], [18])
    @test distance(Count{TransitionMutation}, dnas2) == ([4], [18])
    @test distance(Count{TransversionMutation}, dnas2) == ([8], [18])
    @test distance(Count{Kimura80}, dnas2) == ([4], [8], [18])
    @test distance(Count{AnyMutation}, m2) == ([12], [18])
    @test distance(Count{TransitionMutation}, m2) == ([4], [18])
    @test distance(Count{TransversionMutation}, m2) == ([8], [18])
    @test distance(Count{Kimura80}, m2) == ([4], [8], [18])

    d = distance(Proportion{AnyMutation}, dnas2, 5, 5)
    a = [0.4, 0.8, 1.0, 0.6]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-4
    end
    @test d[2][:] == [5, 5, 3, 5]
    d = distance(Proportion{TransitionMutation}, dnas2, 5, 5)
    a = [0.0, 0.4, 0.333333, 0.2]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-4
    end
    @test d[2][:] == [5, 5, 3, 5]
    d = distance(Proportion{TransversionMutation}, dnas2, 5, 5)
    a = [0.4, 0.4, 0.666667, 0.4]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-4
    end
    @test d[2][:] == [5, 5, 3, 5]

    @test distance(Proportion{AnyMutation}, dnas1) == ([(12 / 16)], [16])
    @test distance(Proportion{TransitionMutation}, dnas1) == ([(4 / 16)], [16])
    @test distance(Proportion{TransversionMutation}, dnas1) == ([(8 / 16)], [16])
    @test distance(Proportion{AnyMutation}, m1) == ([(12 / 16)], [16])
    @test distance(Proportion{TransitionMutation}, m1) == ([(4 / 16)], [16])
    @test distance(Proportion{TransversionMutation}, m1) == ([(8 / 16)], [16])

    @test distance(Proportion{AnyMutation}, dnas2) == ([(12 / 18)], [18])
    @test distance(Proportion{TransitionMutation}, dnas2) == ([(4 / 18)], [18])
    @test distance(Proportion{TransversionMutation}, dnas2) == ([(8 / 18)], [18])
    @test distance(Proportion{AnyMutation}, m2) == ([(12 / 18)], [18])
    @test distance(Proportion{TransitionMutation}, m2) == ([(4 / 18)], [18])
    @test distance(Proportion{TransversionMutation}, m2) == ([(8 / 18)], [18])

    @test distance(JukesCantor69, dnas1) == ([Inf], [Inf]) # Returns infinity as 12/16 is 0.75 - mutation saturation.
    @test distance(JukesCantor69, m1) == ([Inf], [Inf])

    @test round(distance(JukesCantor69, dnas2)[1][1], 3) == 1.648
    @test round(distance(JukesCantor69, dnas2)[2][1], 3) == 1
    @test round(distance(JukesCantor69, m2)[1][1], 3) == 1.648
    @test round(distance(JukesCantor69, m2)[2][1], 3) == 1
    @test_throws DomainError distance(JukesCantor69, dnas2, 5, 5)
    d = distance(JukesCantor69, dnas3, 5, 5)
    a = [0.232616, 0.571605, 0.44084, 0.571605]
    v = [0.0595041, 0.220408, 0.24, 0.220408]
    for i in 1:length(d[1])
        @test_approx_eq_eps d[1][i] a[i] 1e-5
        @test_approx_eq_eps d[2][i] v[i] 1e-5
    end

    @test round(distance(Kimura80, dnas2)[1][1], 3) == 1.648
    @test round(distance(Kimura80, dnas2)[2][1], 3) == 1
    @test round(distance(Kimura80, m2)[1][1], 3) == 1.648
    @test round(distance(Kimura80, m2)[2][1], 3) == 1
end
=#

@testset "VCF" begin
    metainfo = VCF.MetaInfo()
    @test !isfilled(metainfo)
    @test ismatch(r"^Bio.Var.VCF.MetaInfo: <not filled>", repr(metainfo))
    @test_throws ArgumentError metainfotag(metainfo)

    metainfo = VCF.MetaInfo(b"##source=foobar1234")
    @test isfilled(metainfo)
    @test metainfotag(metainfo) == "source"
    @test metainfoval(metainfo) == "foobar1234"

    metainfo = VCF.MetaInfo("##source=foobar1234")
    @test isfilled(metainfo)
    @test metainfotag(metainfo) == "source"
    @test metainfoval(metainfo) == "foobar1234"

    metainfo = VCF.MetaInfo(metainfo)
    @test isa(metainfo, VCF.MetaInfo)
    metainfo = VCF.MetaInfo(metainfo, tag="date")
    @test metainfotag(metainfo) == "date"
    metainfo = VCF.MetaInfo(metainfo, value="2017-01-30")
    @test metainfoval(metainfo) == "2017-01-30"
    metainfo = VCF.MetaInfo(metainfo, tag="INFO", value=["ID"=>"DP", "Number"=>"1", "Type"=>"Integer", "Description"=>"Total Depth"])
    @test metainfo["ID"] == "DP"
    @test metainfo["Number"] == "1"
    @test metainfo["Type"] == "Integer"
    @test metainfo["Description"] == "Total Depth"
    @test metainfotag(metainfo) == "INFO"
    @test metainfoval(metainfo) == """<ID=DP,Number=1,Type=Integer,Description="Total Depth">"""

    record = VCF.Record()
    @test !isfilled(record)
    @test ismatch(r"^Bio.Var.VCF.Record: <not filled>", repr(record))
    @test_throws ArgumentError VCF.chrom(record)

    record = VCF.Record("20\t302\t.\tT\tTA\t999\t.\t.\tGT")
    @test isfilled(record)
    @test VCF.haschrom(record)
    @test VCF.chrom(record) == "20"
    @test VCF.haspos(record)
    @test VCF.pos(record) == 302
    @test !VCF.hasid(record)
    @test_throws MissingFieldException VCF.id(record)
    @test VCF.hasref(record)
    @test VCF.ref(record) == "T"
    @test VCF.hasalt(record)
    @test VCF.alt(record) == ["TA"]
    @test VCF.hasqual(record)
    @test VCF.qual(record) == 999
    @test !VCF.hasfilter(record)
    @test_throws MissingFieldException VCF.filter(record)
    @test VCF.infokeys(record) == String[]
    @test !VCF.hasinfo(record)
    @test_throws MissingFieldException VCF.info(record)
    @test VCF.hasformat(record)
    @test VCF.format(record) == ["GT"]

    # empty data is not a valid VCF record
    @test_throws ArgumentError VCF.Record("")
    @test_throws ArgumentError VCF.Record(b"")

    record = VCF.Record(b".\t.\t.\t.\t.\t.\t.\t.\t")
    @test isfilled(record)
    @test !VCF.haschrom(record)
    @test !VCF.haspos(record)
    @test !VCF.hasid(record)
    @test !VCF.hasref(record)
    @test !VCF.hasalt(record)
    @test !VCF.hasqual(record)
    @test !VCF.hasfilter(record)
    @test !VCF.hasinfo(record)

    record = VCF.Record(record)
    @test isa(record, VCF.Record)
    record = VCF.Record(record, chrom="chr1")
    @test VCF.chrom(record) == "chr1"
    record = VCF.Record(record, pos=1234)
    @test VCF.pos(record) == 1234
    record = VCF.Record(record, id="rs1111")
    @test VCF.id(record) == ["rs1111"]
    record = VCF.Record(record, ref="A")
    @test VCF.ref(record) == "A"
    record = VCF.Record(record, alt=["AT"])
    @test VCF.alt(record) == ["AT"]
    record = VCF.Record(record, qual=11.2)
    @test VCF.qual(record) == 11.2
    record = VCF.Record(record, filter="PASS")
    @test VCF.filter(record) == ["PASS"]
    record = VCF.Record(record, info=Dict("DP" => 20, "AA" => "AT", "DB"=>nothing))
    @test VCF.info(record, "DP") == "20"
    @test VCF.info(record, "AA") == "AT"
    @test VCF.info(record, "DB") == ""
    @test VCF.infokeys(record) == ["DP", "AA", "DB"]
    record = VCF.Record(record, genotype=[Dict("GT" => "0/0", "DP" => [10,20])])
    @test VCF.format(record) == ["DP", "GT"]
    @test VCF.genotype(record) == [["10,20", "0/0"]]

    let header = VCF.Header()
        @test isempty(header)
        push!(header, "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta")
        @test !isempty(header)
        @test length(header) == 1
        unshift!(header, "##fileformat=VCFv4.3")
        @test length(header) == 2
        @test collect(header) == [
            VCF.MetaInfo("##fileformat=VCFv4.3"),
            VCF.MetaInfo("##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta")]
        @test startswith(repr(header), "Bio.Var.VCF.Header:")
    end

    let header = VCF.Header(["##fileformat=VCFv4.3"], ["Sample1"])
        @test !isempty(header)
        @test length(header) == 1
        @test header.sampleID == ["Sample1"]
        @test first(header) == VCF.MetaInfo("##fileformat=VCFv4.3")
    end

    # minimum header
    data = b"""
    ##fileformat=VCFv4.3
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    """
    reader = VCF.Reader(BufferedInputStream(data))
    @test isa(header(reader), VCF.Header)
    let header = header(reader)
        @test length(header.metainfo) == 1
        @test metainfotag(header.metainfo[1]) == "fileformat"
        @test metainfoval(header.metainfo[1]) == "VCFv4.3"
        @test isempty(header.sampleID)
    end

    # realistic header
    data = b"""
    ##fileformat=VCFv4.2
    ##fileDate=20090805
    ##source=myImputationProgramV3.1
    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
    ##phasing=partial
    ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
    """
    reader = VCF.Reader(BufferedInputStream(data))
    @test isa(header(reader), VCF.Header)

    let header = header(reader)
        @test length(header.metainfo) == 10

        let metainfo = header.metainfo[1]
            @test metainfotag(metainfo) == "fileformat"
            @test metainfoval(metainfo) == "VCFv4.2"
            @test_throws ArgumentError keys(metainfo)
            @test_throws ArgumentError values(metainfo)
        end
        @test length(find(header, "fileformat")) == 1
        @test first(find(header, "fileformat")) == header.metainfo[1]

        let metainfo = header.metainfo[2]
            @test metainfotag(metainfo) == "fileDate"
            @test metainfoval(metainfo) == "20090805"
            @test_throws ArgumentError keys(metainfo)
            @test_throws ArgumentError values(metainfo)
        end

        let metainfo = header.metainfo[5]
            @test metainfotag(metainfo) == "contig"
            @test metainfoval(metainfo) == """<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>"""
            @test keys(metainfo) == ["ID", "length", "assembly", "md5", "species", "taxonomy"]
            @test values(metainfo) == ["20", "62435964", "B36", "f126cdf8a6e0c7f379d618ff66beb2da", "Homo sapiens", "x"]
            @test metainfo["ID"] == "20"
            @test metainfo["md5"] == "f126cdf8a6e0c7f379d618ff66beb2da"
            @test metainfo["taxonomy"] == "x"
        end

        let metainfo = header.metainfo[7]
            @test metainfotag(metainfo) == "INFO"
            @test metainfoval(metainfo) == """<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">"""
            @test keys(metainfo) == ["ID", "Number", "Type", "Description"]
            @test values(metainfo) == ["NS", "1", "Integer", "Number of Samples With Data"]
            @test metainfo["ID"] == "NS"
            @test metainfo["Type"] == "Integer"
        end
        @test length(find(header, "INFO")) == 4

        @test header.sampleID == ["NA00001", "NA00002", "NA00003"]
    end

    data = b"""
    ##fileformat=VCFv4.3
    ##contig=<ID=chr1>
    ##contig=<ID=chr2>
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##FORMAT=<ID=GT,Number=1,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002
    chr1\t1234\trs001234\tA\tC\t30\tPASS\tDP=10;AF=0.3\tGT\t0|0\t0/1
    chr2\t4\t.\tA\tAA,AAT\t.\t.\tDP=5\tGT:DP\t0|1:42\t0/1
    """
    reader = VCF.Reader(BufferedInputStream(data))
    record = VCF.Record()

    @test read!(reader, record) === record
    @test VCF.chrom(record) == "chr1"
    @test VCF.pos(record) === 1234
    @test VCF.id(record) == ["rs001234"]
    @test VCF.ref(record) == "A"
    @test VCF.alt(record) == ["C"]
    @test VCF.qual(record) === 30.0
    @test VCF.filter(record) == ["PASS"]
    @test VCF.info(record) == ["DP" => "10", "AF" => "0.3"]
    @test VCF.info(record, "DP") == "10"
    @test VCF.info(record, "AF") == "0.3"
    @test VCF.format(record) == ["GT"]
    @test VCF.genotype(record) == [["0|0"], ["0/1"]]
    @test VCF.genotype(record, 1) == VCF.genotype(record)[1]
    @test VCF.genotype(record, 2) == VCF.genotype(record)[2]
    @test VCF.genotype(record, 1, "GT") == "0|0"
    @test VCF.genotype(record, 2, "GT") == "0/1"
    @test VCF.genotype(record, 1:2, "GT") == ["0|0", "0/1"]
    @test VCF.genotype(record, :, "GT") == VCF.genotype(record, 1:2, "GT")
    @test ismatch(r"^Bio.Var.VCF.Record:\n.*", repr(record))

    @test read!(reader, record) === record
    @test VCF.chrom(record) == "chr2"
    @test VCF.pos(record) == 4
    @test !VCF.hasid(record)
    @test VCF.ref(record) == "A"
    @test VCF.alt(record) == ["AA", "AAT"]
    @test !VCF.hasqual(record)
    @test !VCF.hasfilter(record)
    @test VCF.info(record) == ["DP" => "5"]
    @test VCF.info(record, "DP") == "5"
    @test_throws KeyError VCF.info(record, "AF")
    @test VCF.format(record) == ["GT", "DP"]
    @test VCF.genotype(record) == [["0|1", "42"], ["0/1", "."]]
    @test VCF.genotype(record, 1) == VCF.genotype(record)[1]
    @test VCF.genotype(record, 2) == VCF.genotype(record)[2]
    @test VCF.genotype(record, 1, "GT") == "0|1"
    @test VCF.genotype(record, 1, "DP") == "42"
    @test VCF.genotype(record, 2, "GT") == "0/1"
    @test VCF.genotype(record, 2, "DP") == "."
    @test VCF.genotype(record, 1:2, "GT") == ["0|1", "0/1"]
    @test VCF.genotype(record, 1:2, "DP") == ["42", "."]
    @test VCF.genotype(record, :, "DP") == VCF.genotype(record, 1:2, "DP")
    @test_throws KeyError VCF.genotype(record, :, "BAD")

    @test_throws EOFError read!(reader, record)

    # round-trip test
    vcfdir = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "VCF")
    for specimen in YAML.load_file(joinpath(vcfdir, "index.yml"))
        filepath = joinpath(vcfdir, specimen["filename"])
        records = VCF.Record[]
        reader = open(VCF.Reader, filepath)
        output = IOBuffer()
        writer = VCF.Writer(output, header(reader))
        for record in reader
            write(writer, record)
            push!(records, record)
        end
        close(reader)
        flush(writer)

        records2 = VCF.Record[]
        for record in VCF.Reader(IOBuffer(takebuf_array(output)))
            push!(records2, record)
        end
        @test records == records2
    end
end

function parsehex(str)
    return map(x -> parse(UInt8, x, 16), split(str, ' '))
end

@testset "BCF" begin
    record = BCF.Record()
    @test !isfilled(record)
    @test ismatch(r"^Bio.Var.BCF.Record: <not filled>", repr(record))
    @test_throws ArgumentError BCF.chrom(record)

    record = BCF.Record()
    record.sharedlen = 0x1c
    record.indivlen = 0x00
    # generated from bcftools 1.3.1 (htslib 1.3.1)
    record.data = parsehex("00 00 00 00 ff ff ff ff 01 00 00 00 01 00 80 7f 00 00 01 00 00 00 00 00 07 17 2e 00")
    record.filled = 1:endof(record.data)
    @test BCF.chrom(record) == 1
    record = BCF.Record(record)
    @test isa(record, BCF.Record)
    record = BCF.Record(record, chrom=4)
    @test BCF.chrom(record) == 4
    record = BCF.Record(record, pos=1234)
    @test BCF.pos(record) == 1234
    record = BCF.Record(record, qual=12.3)
    @test BCF.qual(record) == 12.3f0
    record = BCF.Record(record, id="rs1234")
    @test BCF.id(record) == "rs1234"
    record = BCF.Record(record, ref="AT")
    @test BCF.ref(record) == "AT"
    record = BCF.Record(record, alt=["ATT", "ACT"])
    @test BCF.alt(record) == ["ATT", "ACT"]
    record = BCF.Record(record, filter=[2, 3])
    @test BCF.filter(record) == [2, 3]
    record = BCF.Record(record, info=Dict(1 => Int8[42]))
    @test BCF.info(record) == [(1, 42)]
    @test BCF.info(record, simplify=false) == [(1, [42])]
    @test BCF.info(record, 1) == 42
    @test BCF.info(record, 1, simplify=false) == [42]

    bcfdir = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "BCF")
    reader = BCF.Reader(open(joinpath(bcfdir, "example.bcf")))
    let header = header(reader)
        @test length(find(header, "fileformat")) == 1
        @test find(header, "fileformat")[1] == VCF.MetaInfo("##fileformat=VCFv4.2")
        @test length(find(header, "FORMAT")) == 4
    end
    record = BCF.Record()
    @test read!(reader, record) === record
    @test BCF.chrom(record) == 1
    @test BCF.pos(record) == 14370
    @test BCF.rlen(record) == 1
    @test BCF.id(record) == "rs6054257"
    @test BCF.ref(record) == "G"
    @test BCF.alt(record) == ["A"]
    @test BCF.qual(record) == 29.0
    @test BCF.filter(record) == [1]
    @test BCF.info(record) == [(1,3),(2,14),(3,0.5),(5,nothing),(6,nothing)]
    @test BCF.info(record, simplify=false) == [(1,[3]),(2,[14]),(3,[0.5]),(5,[]),(6,[])]
    @test BCF.genotype(record) == [(9,[[2,3],[4,3],[4,4]]),(10,[[48],[48],[43]]),(2,[[1],[8],[5]]),(11,[[51,51],[51,51],[-128,-128]])]
    @test BCF.genotype(record, 1) == [(9, [2,3]), (10, [48]), (2, [1]), (11, [51,51])]
    @test BCF.genotype(record, 1, 9) == [2,3]
    @test BCF.genotype(record, 1, 10) == [48]
    @test BCF.genotype(record, 2, 9) == [4,3]
    @test BCF.genotype(record, :, 9) == [[2,3],[4,3],[4,4]]
    @test ismatch(r"^Bio.Var.BCF.Record:\n.*", repr(record))
    close(reader)

    # round-trip test
    for specimen in YAML.load_file(joinpath(bcfdir, "index.yml"))
        filepath = joinpath(bcfdir, specimen["filename"])
        records = BCF.Record[]
        reader = open(BCF.Reader, filepath)
        output = IOBuffer()
        writer = BCF.Writer(output, header(reader))
        for record in reader
            write(writer, record)
            push!(records, record)
        end
        # HACK: take the data buffer before closing the writer
        data = output.data
        close(reader)
        close(writer)

        records2 = BCF.Record[]
        for record in BCF.Reader(IOBuffer(data))
            push!(records2, record)
        end
        @test records == records2
    end
end

end # module TestVar
