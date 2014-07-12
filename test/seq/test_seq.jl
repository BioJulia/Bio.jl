module TestSeq

using FactCheck
using Bio
using Bio.Seq

# Check that sequences in strings survive the round trip:
#   String → NucleotideSequence → String
function check_string_construction(T::Type, seq::String)
    return convert(String, NucleotideSequence{T}(seq)) == uppercase(seq)
end


# Return a random DNA/RNA sequence of the given length
function random_dna(n::Integer)
    nts = ['A', 'C', 'G', 'T', 'N']
    return convert(String, [nts[i] for i in rand(1:5, n)])
end

function random_rna(n::Integer)
    nts = ['A', 'C', 'G', 'U', 'N']
    return convert(String, [nts[i] for i in rand(1:5, n)])
end


facts("String Construction") do
    # Non-nucleotide characters should throw
    @fact_throws DNASequence("ACCNNCATTTTTTAGATXATAG")
    @fact_throws RNASequence("ACCNNCATTTTTTAGATXATAG")

    # Empty
    @fact check_string_construction(DNANucleotide, "") => true
    @fact check_string_construction(RNANucleotide, "") => true

    reps = 10
    for len in [1, 10, 32, 1000, 10000, 100000]
        @fact all([check_string_construction(DNANucleotide, random_dna(len)) for _ in 1:reps]) => true
        @fact all([check_string_construction(RNANucleotide, random_rna(len)) for _ in 1:reps]) => true
    end
end


facts("Subsequence Construction") do

end


end # TestSeq
