module TestVar

using Base.Test

using Bio.Seq
using Bio.Var

@testset "Counting mutations" begin

    # Create a 20bp test DNA sequence pair containing every possible transition (4),
    # every possible transversion (8), and 2 gapped sites and 2 ambiguous sites.
    # This leaves 4 sites non-mutated/conserved.
    dnas = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
    m1 = seqmatrix(dnas, :seq)

    rnas = [rna"AUUG-ACCUGGNUUUCCGAA", rna"A-ACAGAGUAUACRGUCGUC"]
    m2 = seqmatrix(rnas, :seq)

    @test count_mutations(AnyMutation, dnas) == count_mutations(AnyMutation, rnas) == ([12], [16])
    @test count_mutations(AnyMutation, m1) == count_mutations(AnyMutation, m2) == ([12], [16])
    @test count_mutations(TransitionMutation, dnas) == count_mutations(TransitionMutation, rnas) == ([4], [16])
    @test count_mutations(TransitionMutation, m1) == count_mutations(TransitionMutation, m2) == ([4], [16])
    @test count_mutations(TransversionMutation, dnas) == count_mutations(TransversionMutation, rnas) == ([8], [16])
    @test count_mutations(TransversionMutation, m1) == count_mutations(TransversionMutation, m2) == ([8], [16])
    @test count_mutations(TransitionMutation, TransversionMutation, dnas) == count_mutations(TransitionMutation, TransversionMutation, rnas) == ([4], [8], [16])
    @test count_mutations(TransitionMutation, TransversionMutation, m1) == count_mutations(TransitionMutation, TransversionMutation, m2) == ([4], [8], [16])
    @test count_mutations(TransversionMutation, TransitionMutation, dnas) == count_mutations(TransversionMutation, TransitionMutation, rnas) == ([4], [8], [16])
    @test count_mutations(TransversionMutation, TransitionMutation, m1) == count_mutations(TransversionMutation, TransitionMutation, m2) == ([4], [8], [16])

    ans = Bool[false, false, true, true, false, true, true, true, false, true, true, false, true, false, true, true, false, false, true, true]
    @test flagmutations(AnyMutation, m1)[1][:,1] == ans
    @test flagmutations(AnyMutation, m2)[1][:,1] == ans


end

@testset "Distance Computation" begin

    dnas1 = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
    m1 = seqmatrix(dnas1, :seq)

    dnas2 = [dna"attgaacctggntttccgaa", dna"atacagagtatacrgtcgtc"]
    m2 = seqmatrix(dnas2, :seq)

    @test distance(Count{AnyMutation}, dnas1) == ([12], [16])
    @test distance(Count{TransitionMutation}, dnas1) == ([4], [16])
    @test distance(Count{TransversionMutation}, dnas1) == ([8], [16])
    @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])
    @test distance(Count{AnyMutation}, m1) == ([12], [16])
    @test distance(Count{TransitionMutation}, m1) == ([4], [16])
    @test distance(Count{TransversionMutation}, m1) == ([8], [16])
    @test distance(Count{Kimura80}, m1) == ([4], [8], [16])

    @test distance(Count{AnyMutation}, dnas2) == ([12], [18])
    @test distance(Count{TransitionMutation}, dnas2) == ([4], [18])
    @test distance(Count{TransversionMutation}, dnas2) == ([8], [18])
    @test distance(Count{Kimura80}, dnas2) == ([4], [8], [18])
    @test distance(Count{AnyMutation}, m2) == ([12], [18])
    @test distance(Count{TransitionMutation}, m2) == ([4], [18])
    @test distance(Count{TransversionMutation}, m2) == ([8], [18])
    @test distance(Count{Kimura80}, m2) == ([4], [8], [18])

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

    @test round(distance(Kimura80, dnas2)[1][1], 3) == 1.648
    @test round(distance(Kimura80, dnas2)[2][1], 3) == 1
    @test round(distance(Kimura80, m2)[1][1], 3) == 1.648
    @test round(distance(Kimura80, m2)[2][1], 3) == 1

end

end # module TestVar
