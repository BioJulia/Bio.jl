module TestVar

using Base.Test

using Bio: Seq, Var
using TestFunctions
using PairwiseListMatrices
using IntervalTrees: IntervalValue
import BufferedStreams: BufferedInputStream
import YAML

typealias PWM PairwiseListMatrix

function generate_testcase{A<:Union{DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4}, RNAAlphabet{2}}}(::Type{A}, len::Int)
    a = [convert(Char, i)  for i in alphabet(A)]
    probs = Vector{Float64}(length(a))
    fill!(probs, 1 / length(a))
    return BioSequence{A}(random_seq(len, a, probs))
end

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

@testset "Var" begin
    @testset "Site counting and identification" begin
        @testset "Naive methods" begin

            alphabets = (DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4}, RNAAlphabet{2})

            for alph in alphabets

                # Answers to these tests were worked out manually to verify count_sites_naive was working.
                # seqA and seqB contain all possible observations of sites.

                istwobit = Seq.bitsof(alph) == 2

                seqA, seqB = generate_possibilities_tester(alph)

                # Test methods which work on single sequences.
                @test count_sites_naive(Certain, seqA) == ifelse(istwobit, length(seqA), 49)
                @test count_sites_naive(Certain, seqB) == ifelse(istwobit, length(seqB), 19)
                @test count_sites_naive(Gap, seqA) == ifelse(istwobit, 0, 16)
                @test count_sites_naive(Gap, seqB) == ifelse(istwobit, 0, 1)
                @test count_sites_naive(Ambiguous, seqA) == ifelse(istwobit, 0, length(seqA) - 65)
                @test count_sites_naive(Ambiguous, seqB) == ifelse(istwobit, 0, length(seqB) - 20)

                # Test methods which work on two sequences.
                # Test when sequences are of the same bitencoding.

                @test count_sites_naive(Certain, seqA, seqB) == count_sites_naive(Certain, seqB, seqA) == 10
                @test count_sites_naive(Gap, seqA, seqB) == count_sites_naive(Gap, seqB, seqA) == ifelse(istwobit, 0, 16)
                @test count_sites_naive(Ambiguous, seqA, seqB) == count_sites_naive(Ambiguous, seqB, seqA) == ifelse(istwobit, 0, 121)
                @test count_sites_naive(Match, seqA, seqB) == count_sites_naive(Match, seqB, seqA) == length(alphabet(alph))
                @test count_sites_naive(Mismatch, seqA, seqB) == count_sites_naive(Mismatch, seqB, seqA) == (length(seqA) - length(alphabet(alph)))

                @test count_sites_naive(Mutated, seqA, seqB) == count_sites_naive(Mutated, seqB, seqA) == (6, ifelse(istwobit, 0, 126))
                @test count_sites_naive(Conserved, seqA, seqB) == count_sites_naive(Conserved, seqB, seqA) == (4, ifelse(istwobit, 0, 126))
                @test count_sites_naive(Transition, seqA, seqB) == count_sites_naive(Transition, seqB, seqA) == (2, ifelse(istwobit, 0, 126))
                @test count_sites_naive(Transversion, seqA, seqB) == count_sites_naive(Transversion, seqB, seqA) == (4, ifelse(istwobit, 0, 126))
            end

            # Test for when sequences are of different bitencodings.
            for alphs in [(DNAAlphabet{2}, DNAAlphabet{4}),
                          (RNAAlphabet{2}, RNAAlphabet{4})]
                seqA, seqB = generate_possibilities_tester(alphs...)
                @test count_sites_naive(Certain, seqA, seqB) == count_sites_naive(Certain, seqB, seqA) == 16
                @test count_sites_naive(Gap, seqA, seqB) == count_sites_naive(Gap, seqB, seqA) == 4
                @test count_sites_naive(Ambiguous, seqA, seqB) == count_sites_naive(Ambiguous, seqB, seqA) == 44
                @test count_sites_naive(Match, seqA, seqB) == count_sites_naive(Match, seqB, seqA) == 4
                @test count_sites_naive(Mismatch, seqA, seqB) == count_sites_naive(Mismatch, seqB, seqA) == 60

                @test count_sites_naive(Mutated, seqA, seqB) == count_sites_naive(Mutated, seqB, seqA) == (12, 48)
                @test count_sites_naive(Conserved, seqA, seqB) == count_sites_naive(Conserved, seqB, seqA) == (4, 48)
                @test count_sites_naive(Transition, seqA, seqB) == count_sites_naive(Transition, seqB, seqA) == (4, 48)
                @test count_sites_naive(Transversion, seqA, seqB) == count_sites_naive(Transversion, seqB, seqA) == (8, 48)
            end
        end

        @testset "Pairwise methods" begin
            dnas = [dna"ATCGCCA-", dna"ATCGCCTA", dna"ATCGCCT-", dna"GTCGCCTA"]
            rnas = [rna"AUCGCCA-", rna"AUCGCCUA", rna"AUCGCCU-", rna"GUCGCCUA"]
            answer_mismatch = PWM{Int, false}([0 2 1 3; 2 0 1 1; 1 1 0 2; 3 1 2 0])
            answer_match = PWM{Int, false}([0 6 7 5; 6 0 7 7; 7 7 0 6; 5 7 6 0])
            for i in (dnas, rnas)
                @test count_sites(Mismatch, i) == answer_mismatch
                @test count_sites(Match, i) == answer_match
                @test count_sites(Certain, i) == PWM{Int, false}([0 7 7 7; 7 0 7 8; 7 7 0 7; 7 8 7 0])
                @test count_sites(Ambiguous, i) == PWM{Int, false}([0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0])
                @test count_sites(Gap, i) == PWM{Int, false}([0 1 1 1; 1 0 1 0; 1 1 0 1; 1 0 1 0])
                ambigs = PWM{Int, false}([0 1 1 1; 1 0 1 0; 1 1 0 1; 1 0 1 0])
                @test count_sites(Conserved, i) == (PWM{Int, false}([0 6 6 5; 6 0 7 7; 6 7 0 6; 5 7 6 0]), ambigs)
                @test count_sites(Mutated, i) == (PWM{Int, false}([0 1 1 2; 1 0 0 1; 1 0 0 1; 2 1 1 0]), ambigs)
            end
        end

        @testset "Windowed methods" begin
            dnaA = dna"ATCGCCA-M"
            dnaB = dna"ATCGCCTAA"
            rnaA = rna"AUCGCCA-M"
            rnaB = rna"AUCGCCUAA"
            matches = [3, 3, 3, 3, 2, 1, 0]
            idxes = [1:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9]
            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count_sites(Certain, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                                       IntervalValue(2, 4, 3),
                                                                       IntervalValue(3, 5, 3),
                                                                       IntervalValue(4, 6, 3),
                                                                       IntervalValue(5, 7, 3),
                                                                       IntervalValue(6, 8, 2),
                                                                       IntervalValue(7, 9, 1)]
                @test count_sites(Ambiguous, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                         IntervalValue(2, 4, 0),
                                                                         IntervalValue(3, 5, 0),
                                                                         IntervalValue(4, 6, 0),
                                                                         IntervalValue(5, 7, 0),
                                                                         IntervalValue(6, 8, 0),
                                                                         IntervalValue(7, 9, 1)]
                @test count_sites(Gap, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                   IntervalValue(2, 4, 0),
                                                                   IntervalValue(3, 5, 0),
                                                                   IntervalValue(4, 6, 0),
                                                                   IntervalValue(5, 7, 0),
                                                                   IntervalValue(6, 8, 1),
                                                                   IntervalValue(7, 9, 1)]
                @test count_sites(Match, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                                     IntervalValue(2, 4, 3),
                                                                     IntervalValue(3, 5, 3),
                                                                     IntervalValue(4, 6, 3),
                                                                     IntervalValue(5, 7, 2),
                                                                     IntervalValue(6, 8, 1),
                                                                     IntervalValue(7, 9, 0)]
                @test count_sites(Mismatch, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                        IntervalValue(2, 4, 0),
                                                                        IntervalValue(3, 5, 0),
                                                                        IntervalValue(4, 6, 0),
                                                                        IntervalValue(5, 7, 1),
                                                                        IntervalValue(6, 8, 2),
                                                                        IntervalValue(7, 9, 3)]

                @test count_sites(Conserved, seqs[1], seqs[2], 3, 1) == ([IntervalValue(1, 3, 3),
                                                                          IntervalValue(2, 4, 3),
                                                                          IntervalValue(3, 5, 3),
                                                                          IntervalValue(4, 6, 3),
                                                                          IntervalValue(5, 7, 2),
                                                                          IntervalValue(6, 8, 1),
                                                                          IntervalValue(7, 9, 0)],
                                                                         [IntervalValue(1, 3, 0),
                                                                          IntervalValue(2, 4, 0),
                                                                          IntervalValue(3, 5, 0),
                                                                          IntervalValue(4, 6, 0),
                                                                          IntervalValue(5, 7, 0),
                                                                          IntervalValue(6, 8, 1),
                                                                          IntervalValue(7, 9, 2)])
                @test count_sites(Mutated, seqs[1], seqs[2], 3, 1) == ([IntervalValue(1, 3, 0),
                                                                        IntervalValue(2, 4, 0),
                                                                        IntervalValue(3, 5, 0),
                                                                        IntervalValue(4, 6, 0),
                                                                        IntervalValue(5, 7, 1),
                                                                        IntervalValue(6, 8, 1),
                                                                        IntervalValue(7, 9, 1)],
                                                                       [IntervalValue(1, 3, 0),
                                                                        IntervalValue(2, 4, 0),
                                                                        IntervalValue(3, 5, 0),
                                                                        IntervalValue(4, 6, 0),
                                                                        IntervalValue(5, 7, 0),
                                                                        IntervalValue(6, 8, 1),
                                                                        IntervalValue(7, 9, 2)])
                @test count_sites(Transition, seqs[1], seqs[2], 3, 1) == ([IntervalValue(1, 3, 0),
                                                                           IntervalValue(2, 4, 0),
                                                                           IntervalValue(3, 5, 0),
                                                                           IntervalValue(4, 6, 0),
                                                                           IntervalValue(5, 7, 0),
                                                                           IntervalValue(6, 8, 0),
                                                                           IntervalValue(7, 9, 0)],
                                                                          [IntervalValue(1, 3, 0),
                                                                           IntervalValue(2, 4, 0),
                                                                           IntervalValue(3, 5, 0),
                                                                           IntervalValue(4, 6, 0),
                                                                           IntervalValue(5, 7, 0),
                                                                           IntervalValue(6, 8, 1),
                                                                           IntervalValue(7, 9, 2)])
                @test count_sites(Transversion, seqs[1], seqs[2], 3, 1) == ([IntervalValue(1, 3, 0),
                                                                           IntervalValue(2, 4, 0),
                                                                           IntervalValue(3, 5, 0),
                                                                           IntervalValue(4, 6, 0),
                                                                           IntervalValue(5, 7, 1),
                                                                           IntervalValue(6, 8, 1),
                                                                           IntervalValue(7, 9, 1)],
                                                                          [IntervalValue(1, 3, 0),
                                                                           IntervalValue(2, 4, 0),
                                                                           IntervalValue(3, 5, 0),
                                                                           IntervalValue(4, 6, 0),
                                                                           IntervalValue(5, 7, 0),
                                                                           IntervalValue(6, 8, 1),
                                                                           IntervalValue(7, 9, 2)])
            end
        end

    end
end


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
    metainfo = VCFMetaInfo()
    @test !isfilled(metainfo)
    @test ismatch(r"^Bio.Var.VCFMetaInfo: <not filled>", repr(metainfo))
    @test_throws ArgumentError metainfotag(metainfo)

    metainfo = VCFMetaInfo(b"##source=foobar1234")
    @test isfilled(metainfo)
    @test metainfotag(metainfo) == "source"
    @test metainfoval(metainfo) == "foobar1234"

    metainfo = VCFMetaInfo("##source=foobar1234")
    @test isfilled(metainfo)
    @test metainfotag(metainfo) == "source"
    @test metainfoval(metainfo) == "foobar1234"

    metainfo = VCFMetaInfo(metainfo)
    @test isa(metainfo, VCFMetaInfo)
    metainfo = VCFMetaInfo(metainfo, tag="date")
    @test metainfotag(metainfo) == "date"
    metainfo = VCFMetaInfo(metainfo, value="2017-01-30")
    @test metainfoval(metainfo) == "2017-01-30"
    metainfo = VCFMetaInfo(metainfo, tag="INFO", value=["ID"=>"DP", "Number"=>"1", "Type"=>"Integer", "Description"=>"Total Depth"])
    @test metainfo["ID"] == "DP"
    @test metainfo["Number"] == "1"
    @test metainfo["Type"] == "Integer"
    @test metainfo["Description"] == "Total Depth"
    @test metainfotag(metainfo) == "INFO"
    @test metainfoval(metainfo) == """<ID=DP,Number=1,Type=Integer,Description="Total Depth">"""

    record = VCFRecord()
    @test !isfilled(record)
    @test ismatch(r"^Bio.Var.VCFRecord: <not filled>", repr(record))
    @test_throws ArgumentError chromosome(record)

    record = VCFRecord("20\t302\t.\tT\tTA\t999\t.\t.\tGT")
    @test isfilled(record)
    @test haschromosome(record)
    @test chromosome(record) == "20"
    @test hasleftposition(record)
    @test leftposition(record) == 302
    @test !hasidentifier(record)
    @test_throws MissingFieldException identifier(record)
    @test hasreference(record)
    @test reference(record) == "T"
    @test hasalternate(record)
    @test alternate(record) == ["TA"]
    @test hasquality(record)
    @test quality(record) == 999
    @test !hasfilter_(record)
    @test_throws MissingFieldException filter_(record)
    @test infokeys(record) == String[]
    @test !hasinformation(record)
    @test_throws MissingFieldException information(record)
    @test hasformat(record)
    @test format(record) == ["GT"]

    # empty data is not a valid VCF record
    @test_throws ArgumentError VCFRecord("")
    @test_throws ArgumentError VCFRecord(b"")

    record = VCFRecord(b".\t.\t.\t.\t.\t.\t.\t.\t")
    @test isfilled(record)
    @test !haschromosome(record)
    @test !hasleftposition(record)
    @test !hasidentifier(record)
    @test !hasreference(record)
    @test !hasalternate(record)
    @test !hasquality(record)
    @test !hasfilter_(record)
    @test !hasinformation(record)

    record = VCFRecord(record)
    @test isa(record, VCFRecord)
    record = VCFRecord(record, chromosome="chr1")
    @test chromosome(record) == "chr1"
    record = VCFRecord(record, position=1234)
    @test leftposition(record) == 1234
    record = VCFRecord(record, identifier="rs1111")
    @test identifier(record) == ["rs1111"]
    record = VCFRecord(record, reference="A")
    @test reference(record) == "A"
    record = VCFRecord(record, alternate=["AT"])
    @test alternate(record) == ["AT"]
    record = VCFRecord(record, quality=11.2)
    @test quality(record) == 11.2
    record = VCFRecord(record, filter="PASS")
    @test filter_(record) == ["PASS"]
    record = VCFRecord(record, information=Dict("DP" => 20, "AA" => "AT", "DB"=>nothing))
    @test information(record, "DP") == "20"
    @test information(record, "AA") == "AT"
    @test information(record, "DB") == ""
    @test infokeys(record) == ["DP", "AA", "DB"]
    record = VCFRecord(record, genotype=[Dict("GT" => "0/0", "DP" => [10,20])])
    @test format(record) == ["DP", "GT"]
    @test genotype(record) == [["10,20", "0/0"]]

    let header = VCFHeader()
        @test isempty(header)
        push!(header, "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta")
        @test !isempty(header)
        @test length(header) == 1
        unshift!(header, "##fileformat=VCFv4.3")
        @test length(header) == 2
        @test collect(header) == [
            VCFMetaInfo("##fileformat=VCFv4.3"),
            VCFMetaInfo("##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta")]
        @test startswith(repr(header), "Bio.Var.VCFHeader:")
    end

    let header = VCFHeader(["##fileformat=VCFv4.3"], ["Sample1"])
        @test !isempty(header)
        @test length(header) == 1
        @test header.sampleID == ["Sample1"]
        @test first(header) == VCFMetaInfo("##fileformat=VCFv4.3")
    end

    # minimum header
    data = b"""
    ##fileformat=VCFv4.3
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    """
    reader = VCFReader(BufferedInputStream(data))
    @test isa(header(reader), VCFHeader)
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
    reader = VCFReader(BufferedInputStream(data))
    @test isa(header(reader), VCFHeader)

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
    reader = VCFReader(BufferedInputStream(data))
    record = VCFRecord()

    @test read!(reader, record) === record
    @test chromosome(record) == "chr1"
    @test leftposition(record) === 1234
    @test identifier(record) == ["rs001234"]
    @test reference(record) == "A"
    @test alternate(record) == ["C"]
    @test quality(record) === 30.0
    @test filter_(record) == ["PASS"]
    @test information(record) == ["DP" => "10", "AF" => "0.3"]
    @test information(record, "DP") == "10"
    @test information(record, "AF") == "0.3"
    @test format(record) == ["GT"]
    @test genotype(record) == [["0|0"], ["0/1"]]
    @test genotype(record, 1) == genotype(record)[1]
    @test genotype(record, 2) == genotype(record)[2]
    @test genotype(record, 1, "GT") == "0|0"
    @test genotype(record, 2, "GT") == "0/1"
    @test genotype(record, 1:2, "GT") == ["0|0", "0/1"]
    @test genotype(record, :, "GT") == genotype(record, 1:2, "GT")
    @test ismatch(r"^Bio.Var.VCFRecord:\n.*", repr(record))

    @test read!(reader, record) === record
    @test chromosome(record) == "chr2"
    @test leftposition(record) == 4
    @test !hasidentifier(record)
    @test reference(record) == "A"
    @test alternate(record) == ["AA", "AAT"]
    @test !hasquality(record)
    @test !hasfilter_(record)
    @test information(record) == ["DP" => "5"]
    @test information(record, "DP") == "5"
    @test_throws KeyError information(record, "AF")
    @test format(record) == ["GT", "DP"]
    @test genotype(record) == [["0|1", "42"], ["0/1", "."]]
    @test genotype(record, 1) == genotype(record)[1]
    @test genotype(record, 2) == genotype(record)[2]
    @test genotype(record, 1, "GT") == "0|1"
    @test genotype(record, 1, "DP") == "42"
    @test genotype(record, 2, "GT") == "0/1"
    @test genotype(record, 2, "DP") == "."
    @test genotype(record, 1:2, "GT") == ["0|1", "0/1"]
    @test genotype(record, 1:2, "DP") == ["42", "."]
    @test genotype(record, :, "DP") == genotype(record, 1:2, "DP")
    @test_throws KeyError genotype(record, :, "BAD")

    @test_throws EOFError read!(reader, record)

    # round-trip test
    vcfdir = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "VCF")
    for specimen in YAML.load_file(joinpath(vcfdir, "index.yml"))
        filepath = joinpath(vcfdir, specimen["filename"])
        records = VCFRecord[]
        reader = open(VCFReader, filepath)
        output = IOBuffer()
        writer = VCFWriter(output, header(reader))
        for record in reader
            write(writer, record)
            push!(records, record)
        end
        close(reader)
        flush(writer)

        records2 = VCFRecord[]
        for record in VCFReader(IOBuffer(takebuf_array(output)))
            push!(records2, record)
        end
        @test records == records2
    end
end

function parsehex(str)
    return map(x -> parse(UInt8, x, 16), split(str, ' '))
end

@testset "BCF" begin
    record = BCFRecord()
    @test !isfilled(record)
    @test ismatch(r"^Bio.Var.BCFRecord: <not filled>", repr(record))
    @test_throws ArgumentError chromosome(record)

    record = BCFRecord()
    record.sharedlen = 0x1c
    record.indivlen = 0x00
    # generated from bcftools 1.3.1 (htslib 1.3.1)
    record.data = parsehex("00 00 00 00 ff ff ff ff 01 00 00 00 01 00 80 7f 00 00 01 00 00 00 00 00 07 17 2e 00")
    record.filled = 1:endof(record.data)
    @test chromosome(record) == 1
    record = BCFRecord(record)
    @test isa(record, BCFRecord)
    record = BCFRecord(record, chromosome=4)
    @test chromosome(record) == 4
    record = BCFRecord(record, position=1234)
    @test leftposition(record) == 1234
    record = BCFRecord(record, quality=12.3)
    @test quality(record) == 12.3f0
    record = BCFRecord(record, identifier="rs1234")
    @test identifier(record) == "rs1234"
    record = BCFRecord(record, reference="AT")
    @test reference(record) == "AT"
    record = BCFRecord(record, alternate=["ATT", "ACT"])
    @test alternate(record) == ["ATT", "ACT"]
    record = BCFRecord(record, filter=[2, 3])
    @test filter_(record) == [2, 3]
    record = BCFRecord(record, information=Dict(1 => Int8[42]))
    @test information(record) == [(1, 42)]
    @test information(record, simplify=false) == [(1, [42])]
    @test information(record, 1) == 42
    @test information(record, 1, simplify=false) == [42]

    bcfdir = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "BCF")
    reader = BCFReader(open(joinpath(bcfdir, "example.bcf")))
    let header = header(reader)
        @test length(find(header, "fileformat")) == 1
        @test find(header, "fileformat")[1] == VCFMetaInfo("##fileformat=VCFv4.2")
        @test length(find(header, "FORMAT")) == 4
    end
    record = BCFRecord()
    @test read!(reader, record) === record
    @test chromosome(record) == 1
    @test leftposition(record) == 14370
    @test identifier(record) == "rs6054257"
    @test reference(record) == "G"
    @test alternate(record) == ["A"]
    @test quality(record) == 29.0
    @test filter_(record) == [1]
    @test information(record) == [(1,3),(2,14),(3,0.5),(5,nothing),(6,nothing)]
    @test information(record, simplify=false) == [(1,[3]),(2,[14]),(3,[0.5]),(5,[]),(6,[])]
    @test genotype(record) == [(9,[[2,3],[4,3],[4,4]]),(10,[[48],[48],[43]]),(2,[[1],[8],[5]]),(11,[[51,51],[51,51],[-128,-128]])]
    @test genotype(record, 1) == [(9, [2,3]), (10, [48]), (2, [1]), (11, [51,51])]
    @test genotype(record, 1, 9) == [2,3]
    @test genotype(record, 1, 10) == [48]
    @test genotype(record, 2, 9) == [4,3]
    @test genotype(record, :, 9) == [[2,3],[4,3],[4,4]]
    @test ismatch(r"^Bio.Var.BCFRecord:\n.*", repr(record))
    close(reader)

    # round-trip test
    for specimen in YAML.load_file(joinpath(bcfdir, "index.yml"))
        filepath = joinpath(bcfdir, specimen["filename"])
        records = BCFRecord[]
        reader = open(BCFReader, filepath)
        output = IOBuffer()
        writer = BCFWriter(output, header(reader))
        for record in reader
            write(writer, record)
            push!(records, record)
        end
        # HACK: take the data buffer before closing the writer
        data = output.data
        close(reader)
        close(writer)

        records2 = BCFRecord[]
        for record in BCFReader(IOBuffer(data))
            push!(records2, record)
        end
        @test records == records2
    end
end

end # module TestVar
