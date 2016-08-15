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
dna1 = dna"ATTG-ACCTGGNTTTCCGAA"
dna2 = dna"A-ACAGAGTATACRGTCGTC"

rna1 = rna"AUUG-ACCUGGNUUUCCGAA"
rna2 = rna"A-ACAGAGUAUACRGUCGUC"

@test count_mutations(dna1, dna2, DifferentMutation) == count_mutations(rna1, rna2, DifferentMutation) == (12, 16)
@test count_mutations(dna1, dna2, TransitionMutation) == count_mutations(rna1, rna2, TransitionMutation) == (4, 16)
@test count_mutations(dna1, dna2, TransversionMutation) == count_mutations(rna1, rna2, TransversionMutation) == (8, 16)
@test count_mutations(dna1, dna2, TransitionMutation, TransversionMutation) == count_mutations(rna1, rna2, TransitionMutation, TransversionMutation) == (4, 8, 16)
@test count_mutations(dna1, dna2, TransversionMutation, TransitionMutation) == count_mutations(rna1, rna2, TransversionMutation, TransitionMutation) == (4, 8, 16)

end

end # module TestVar
