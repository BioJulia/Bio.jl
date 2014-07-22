module TestSeq

using FactCheck
using Bio
using Bio.Seq

# Return a random DNA/RNA sequence of the given length
function random_seq(n::Integer, nts, probs)
    cumprobs = cumsum(probs)
    x = Array(Char, n)
    for i in 1:n
        x[i] = nts[searchsorted(cumprobs, rand()).start]
    end
    convert(String, x)
end


function random_dna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'T', 'N'], probs)
end


function random_rna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'U', 'N'], probs)
end


function random_dna_kmer(len)
    return random_dna(len, [0.25, 0.25, 0.25, 0.25])
end


function random_rna_kmer(len)
    return random_rna(len, [0.25, 0.25, 0.25, 0.25])
end



function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
end


facts("Construction") do
    # Non-nucleotide characters should throw
    @fact_throws DNASequence("ACCNNCATTTTTTAGATXATAG")
    @fact_throws RNASequence("ACCNNCATTTTTTAGATXATAG")

    # Check that sequences in strings survive round trip conversion:
    #   String → NucleotideSequence → String
    function check_string_construction(T::Type, seq::String)
        return convert(String, NucleotideSequence{T}(seq)) == uppercase(seq)
    end

    reps = 10
    for len in [0, 1, 10, 32, 1000, 10000, 100000]
        @fact all([check_string_construction(DNANucleotide, random_dna(len)) for _ in 1:reps]) => true
        @fact all([check_string_construction(RNANucleotide, random_rna(len)) for _ in 1:reps]) => true
        @fact all([check_string_construction(DNANucleotide, lowercase(random_dna(len))) for _ in 1:reps]) => true
        @fact all([check_string_construction(RNANucleotide, lowercase(random_rna(len))) for _ in 1:reps]) => true
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

    context("Kmer Construction") do
        # Check that kmers in strings survive round trip conversion:
        #   String → Kmer → String
        function check_string_construction(T::Type, seq::String)
            return convert(String, convert(Kmer{T}, seq)) == uppercase(seq)
        end

        # Check that kmers in strings survive round trip conversion:
        #   String → NucleotideSequence → Kmer → NucleotideSequence → String
        function check_roundabout_construction(T::Type, seq::String)
            return convert(String,
                convert(NucleotideSequence{T},
                    convert(Kmer,
                        convert(NucleotideSequence{T}, seq)))) == uppercase(seq)
        end

        for len in [0, 1, 16, 32]
            check_string_construction(DNANucleotide, random_dna_kmer(len))
            check_string_construction(RNANucleotide, random_rna_kmer(len))
            check_roundabout_construction(DNANucleotide, random_dna_kmer(len))
            check_roundabout_construction(RNANucleotide, random_rna_kmer(len))
        end

        # N is not allowed in Kmers
        @fact_throws dnakmer("ACGTNACGT")
        @fact_throws rnakmer("ACGUNACGU")
    end
end


facts("Transforms") do
    context("Reversal") do
        function check_reversal(T, seq)
            return reverse(seq) == convert(String, reverse(convert(T, seq)))
        end

        reps = 10
        for len in [0, 1, 10, 32, 1000, 10000, 100000]
            @fact all([check_reversal(DNASequence, random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_reversal(RNASequence, random_rna(len)) for _ in 1:reps]) => true
        end

        for len in [0, 1, 16, 32]
            @fact all([check_reversal(DNAKmer, random_dna_kmer(len)) for _ in 1:reps]) => true
            @fact all([check_reversal(RNAKmer, random_rna_kmer(len)) for _ in 1:reps]) => true
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
        function check_dna_complement(T, seq)
            return dna_complement(seq) ==
                convert(String, complement(convert(T, seq)))
        end

        function check_rna_complement(T, seq)
            return rna_complement(seq) ==
                convert(String, complement(convert(T, seq)))
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_dna_complement(DNASequence, random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_rna_complement(RNASequence, random_rna(len)) for _ in 1:reps]) => true
        end

        for len in [0, 1, 16, 32]
            @fact all([check_dna_complement(DNAKmer, random_dna_kmer(len)) for _ in 1:reps]) => true
            @fact all([check_rna_complement(RNAKmer, random_rna_kmer(len)) for _ in 1:reps]) => true
        end
    end

    context("Reverse Complement") do
        function check_dna_revcomp(T, seq)
            return reverse(dna_complement(seq)) ==
                convert(String, reverse_complement(convert(T, seq)))
        end

        function check_rna_revcomp(T, seq)
            return reverse(rna_complement(seq)) ==
                convert(String, reverse_complement(convert(T, seq)))
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_dna_revcomp(DNASequence, random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_rna_revcomp(RNASequence, random_rna(len)) for _ in 1:reps]) => true
        end

        for len in [0, 1, 16, 32]
            @fact all([check_dna_revcomp(DNAKmer, random_dna_kmer(len)) for _ in 1:reps]) => true
            @fact all([check_rna_revcomp(RNAKmer, random_rna_kmer(len)) for _ in 1:reps]) => true
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
            return mismatches(convert(T, a), convert(T, b)) == count
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_mismatches(DNASequence, random_dna(len), random_dna(len))
                       for _ in 1:reps]) => true
            @fact all([check_mismatches(RNASequence, random_rna(len), random_rna(len))
                       for _ in 1:reps]) => true
        end

        for len in [0, 1, 16, 32]
            @fact all([check_mismatches(DNAKmer, random_dna_kmer(len), random_dna_kmer(len))
                      for _ in 1:reps]) => true
            @fact all([check_mismatches(RNAKmer, random_rna_kmer(len), random_rna_kmer(len))
                      for _ in 1:reps]) => true
        end
    end
end


facts("Iteration") do
    context("NPositions") do
        function check_ns(T, seq)
            expected = Int[]
            for i in 1:length(seq)
                if seq[i] == 'N'
                    push!(expected, i)
                end
            end

            collect(npositions(T(seq))) == expected
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_ns(DNASequence, random_dna(len)) for _ in 1:reps]) => true
            @fact all([check_ns(RNASequence, random_rna(len)) for _ in 1:reps]) => true
        end
    end

    context("EachKmer") do
        function string_eachkmer(seq::String, k)
            kmers = String[]
            for i in 1:(length(seq) - k + 1)
                subseq = seq[i:i + k - 1]
                if !in('N', subseq)
                    push!(kmers, subseq)
                end
            end
            return kmers
        end

        function check_eachkmer(T, seq::String, k)
            xs = [convert(String, x) for x in collect(eachkmer(NucleotideSequence{T}(seq), k))]
            ys = string_eachkmer(seq, k)
            return xs == ys
        end

        reps = 10
        len = 10000
        for k in [1, 16, 32]
            @fact all([check_eachkmer(DNANucleotide, random_dna(len), k) for _ in 1:reps]) => true
            @fact all([check_eachkmer(RNANucleotide, random_rna(len), k) for _ in 1:reps]) => true
        end

        @fact isempty(collect(eachkmer(dna"", 1))) => true
        @fact isempty(collect(eachkmer(dna"NNNNNNNNNN", 1))) => true
        @fact isempty(collect(eachkmer(dna"ACGT", 0))) => true
        @fact_throws eachkmer(dna"ACGT", -1)
        @fact_throws eachkmer(dna"ACGT", 33)
    end
end


end # TestSeq
