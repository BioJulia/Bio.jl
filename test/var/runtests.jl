module TestVar

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Seq
using Bio.Var

@testset "Counting mutations" begin

    # Create a 20bp test DNA sequence pair containing every possible transition (4),
    # every possible transversion (8), and 2 gapped sites and 2 ambiguous sites.
    # This leaves 4 sites non-mutated/conserved.
    dnas = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]

    rnas = [rna"AUUG-ACCUGGNUUUCCGAA", rna"A-ACAGAGUAUACRGUCGUC"]

    @test count_mutations(dnas, DifferentMutation) == count_mutations(rnas, DifferentMutation) == ([12], [16])
    @test count_mutations(dnas, TransitionMutation) == count_mutations(rnas, TransitionMutation) == ([4], [16])
    @test count_mutations(dnas, TransversionMutation) == count_mutations(rnas, TransversionMutation) == ([8], [16])
    @test count_mutations(dnas, TransitionMutation, TransversionMutation) == count_mutations(rnas, TransitionMutation, TransversionMutation) == ([4], [8], [16])
    @test count_mutations(dnas, TransversionMutation, TransitionMutation) == count_mutations(rnas, TransversionMutation, TransitionMutation) == ([4], [8], [16])

end

# @testset "Distance Computation" begin
#
#     dna1 = dna"ATTG-ACCTGGNTTTCCGAA"
#     dna2 = dna"A-ACAGAGTATACRGTCGTC"
#
#     dna3 = dna"attgaacctggntttccgaa"
#     dna4 = dna"atacagagtatacrgtcgtc"
#
#     @test distance(Count{DifferentMutation}, dna1, dna2) == (12, 16)
#     @test distance(Count{TransitionMutation}, dna1, dna2) == (4, 16)
#     @test distance(Count{TransversionMutation}, dna1, dna2) == (8, 16)
#     @test distance(Count{Kimura80}, dna1, dna2) == (4, 8, 16)
#
#     @test distance(Count{DifferentMutation}, dna3, dna4) == (12, 18)
#     @test distance(Count{TransitionMutation}, dna3, dna4) == (4, 18)
#     @test distance(Count{TransversionMutation}, dna3, dna4) == (8, 18)
#     @test distance(Count{Kimura80}, dna3, dna4) == (4, 8, 18)
#
#     @test distance(Proportion{DifferentMutation}, dna1, dna2) == ((12 / 16), 16)
#     @test distance(Proportion{TransitionMutation}, dna1, dna2) == ((4 / 16), 16)
#     @test distance(Proportion{TransversionMutation}, dna1, dna2) == ((8 / 16), 16)
#
#     @test distance(Proportion{DifferentMutation}, dna3, dna4) == ((12 / 18), 18)
#     @test distance(Proportion{TransitionMutation}, dna3, dna4) == ((4 / 18), 18)
#     @test distance(Proportion{TransversionMutation}, dna3, dna4) == ((8 / 18), 18)
#
#     @test distance(JukesCantor69, dna1, dna2) == (Inf, Inf) # Returns infinity as 12/16 is 0.75 - mutation saturation.
#
#     @test round(distance(JukesCantor69, dna3, dna4)[1], 3) == 1.648
#     @test round(distance(JukesCantor69, dna3, dna4)[2], 3) == 1
#
#     @test round(distance(Kimura80, dna3, dna4)[1], 3) == 1.648
#     @test round(distance(Kimura80, dna3, dna4)[2], 3) == 1
#
# end

end # module TestVar
