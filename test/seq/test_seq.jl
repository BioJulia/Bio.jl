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

function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
end


facts("Construction") do
    # Non-nucleotide characters should throw
    @fact_throws DNASequence("ACCNNCATTTTTTAGATXATAG")
    @fact_throws RNASequence("ACCNNCATTTTTTAGATXATAG")

    reps = 10
    for len in [0, 1, 10, 32, 1000, 10000, 100000]
        @fact all([check_string_construction(DNANucleotide, random_dna(len)) for _ in 1:reps]) => true
        @fact all([check_string_construction(RNANucleotide, random_rna(len)) for _ in 1:reps]) => true
    end

    context("Copy") do
        function check_copy(T, seq)
            return convert(String, copy(NucleotideSequence{T}(seq))) == seq
        end

        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_copy(DNANucleotide, random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_copy(RNANucleotide, random_rna(len)) for _ in 1:reps]) => true
        end
    end

    context("Subsequence Construction") do
        for len in [1, 10, 32, 1000, 10000, 100000]
            seq = random_dna(len)
            dnaseq = DNASequence(seq)

            results = Bool[]
            for _ in 1:reps
                part = random_interval(1, length(seq))
                push!(results, seq[part] == convert(String, dnaseq[part]))
            end
            @fact all(results) => true
        end

        for len in [1, 10, 32, 1000, 10000, 100000]
            seq = random_rna(len)
            rnaseq = RNASequence(seq)

            results = Bool[]
            for _ in 1:reps
                part = random_interval(1, length(seq))

                push!(results, seq[part] == convert(String, rnaseq[part]))
            end
            @fact all(results) => true
        end
    end
end


facts("Transforms") do
    context("Reversal") do
        function check_reversal(T, seq)
            return reverse(seq) == convert(String, reverse(T(seq)))
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_reversal(DNASequence, random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_reversal(RNASequence, random_rna(len)) for _ in 1:reps]) => true
        end
    end

    function dna_complement(seq::String)
        seqc = Array(Char, length(seq))
        for (i, c) in enumerate(seq)
            if c == 'A'
                seqc[i] = 'T'
            elseif c == 'C'
                seqc[i] = 'G'
            elseif c == 'G'
                seqc[i] = 'C'
            elseif c == 'T'
                seqc[i] = 'A'
            else
                seqc[i] = 'N'
            end
        end

        return convert(String, seqc)
    end

    function rna_complement(seq::String)
        seqc = Array(Char, length(seq))
        for (i, c) in enumerate(seq)
            if c == 'A'
                seqc[i] = 'U'
            elseif c == 'C'
                seqc[i] = 'G'
            elseif c == 'G'
                seqc[i] = 'C'
            elseif c == 'U'
                seqc[i] = 'A'
            else
                seqc[i] = 'N'
            end
        end

        return convert(String, seqc)
    end

    context("Complement") do
        function check_dna_complement(seq)
            return dna_complement(seq) == convert(String, complement(DNASequence(seq)))
        end

        function check_rna_complement(seq)
            return rna_complement(seq) == convert(String, complement(RNASequence(seq)))
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_dna_complement(random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_rna_complement(random_rna(len)) for _ in 1:reps]) => true
        end
    end

    context("Reverse Complement") do
        function check_dna_revcomp(seq)
            return reverse(dna_complement(seq)) ==
                convert(String, reverse_complement(DNASequence(seq)))
        end

        function check_rna_revcomp(seq)
            return reverse(rna_complement(seq)) ==
                convert(String, reverse_complement(RNASequence(seq)))
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_dna_revcomp(random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_rna_revcomp(random_rna(len)) for _ in 1:reps]) => true
        end
    end
end


facts("Compare") do
    context("Mismatches") do
        function check_mismatches(T, a, b)
            count = 0
            for (ca, cb) in zip(a, b)
                if ca != cb
                    count += 1
                end
            end
            if mismatches(T(a), T(b)) != count
                println(STDERR, a, "\n", b, "\n", mismatches(T(a), T(b)), "\n", count)
            end

            return mismatches(T(a), T(b)) == count
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_mismatches(DNASequence, random_dna(len), random_dna(len))
                       for _ in 1:reps]) => true
            @fact all([check_mismatches(RNASequence, random_rna(len), random_rna(len))
                       for _ in 1:reps]) => true
        end
    end
end


facts("Iteration") do
    context("Ns") do
        function check_ns(T, seq)
            expected = Int[]
            for i in 1:length(seq)
                if seq[i] == 'N'
                    push!(expected, i)
                end
            end

            collect(ns(T(seq))) == expected
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_ns(DNASequence, random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_ns(RNASequence, random_rna(len)) for _ in 1:reps]) => true
        end
    end
end



end # TestSeq
