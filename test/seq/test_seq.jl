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

const codons = [
        "AAA", "AAC", "AAG", "AAU",
        "ACA", "ACC", "ACG", "ACU",
        "AGA", "AGC", "AGG", "AGU",
        "AUA", "AUC", "AUG", "AUU",
        "CAA", "CAC", "CAG", "CAU",
        "CCA", "CCC", "CCG", "CCU",
        "CGA", "CGC", "CGG", "CGU",
        "CUA", "CUC", "CUG", "CUU",
        "GAA", "GAC", "GAG", "GAU",
        "GCA", "GCC", "GCG", "GCU",
        "GGA", "GGC", "GGG", "GGU",
        "GUA", "GUC", "GUG", "GUU",
        "UAC", "UAU", "UCA", "UCC",
        "UCG", "UCU", "UGC", "UGG",
        "UGU", "UUA", "UUC", "UUG",
        "UUU" ]


function random_translatable_rna(n)
    probs = fill(1.0 / length(codons), length(codons))
    cumprobs = cumsum(probs)
    r = rand()
    x = Array(String, n)
    for i in 1:n
        x[i] = codons[searchsorted(cumprobs, rand()).start]
    end

    return string(x...)
end


function random_dna_kmer(len)
    return random_dna(len, [0.25, 0.25, 0.25, 0.25])
end


function random_rna_kmer(len)
    return random_rna(len, [0.25, 0.25, 0.25, 0.25])
end


function random_aa(len)
    return random_seq(len,
        ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X' ],
        push!(fill(0.049, 20), 0.02))
end


function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
end


facts("NucleotideSequence Construction") do
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


facts("AminoAcidSequence Construction") do
    # Non-aa characters should throw
    @fact_throws AminoAcidSequence("ATGHLMYZZACAGNM")

    # Check that sequences in strings survive round trip conversion:
    #   String → AminoAcidSequence → String
    function check_string_construction(seq::String)
        return convert(String, AminoAcidSequence(seq)) == uppercase(seq)
    end

    reps = 10
    for len in [0, 1, 10, 32, 1000, 10000, 100000]
        @fact all([check_string_construction(random_aa(len)) for _ in 1:reps]) => true
        @fact all([check_string_construction(lowercase(random_aa(len))) for _ in 1:reps]) => true
    end

    context("Copy") do
        function check_copy(seq)
            return convert(String, copy(AminoAcidSequence(seq))) == uppercase(seq)
        end

        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_copy(random_aa(len)) for _ in 1:reps]) => true
        end
    end

    context("Subsequence Construction") do
        for len in [1, 10, 32, 1000, 10000, 100000]
            seq = random_aa(len)
            aaseq = AminoAcidSequence(seq)

            results = Bool[]
            for _ in 1:reps
                part = random_interval(1, length(seq))
                push!(results, seq[part] == convert(String, aaseq[part]))
            end
            @fact all(results) => true
        end
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

    context("Translate") do
        # crummy string translation to test against
        standard_genetic_code_dict = [
            "AAA" => 'K', "AAC" => 'N', "AAG" => 'K', "AAU" => 'N',
            "ACA" => 'T', "ACC" => 'T', "ACG" => 'T', "ACU" => 'T',
            "AGA" => 'R', "AGC" => 'S', "AGG" => 'R', "AGU" => 'S',
            "AUA" => 'I', "AUC" => 'I', "AUG" => 'M', "AUU" => 'I',
            "CAA" => 'Q', "CAC" => 'H', "CAG" => 'Q', "CAU" => 'H',
            "CCA" => 'P', "CCC" => 'P', "CCG" => 'P', "CCU" => 'P',
            "CGA" => 'R', "CGC" => 'R', "CGG" => 'R', "CGU" => 'R',
            "CUA" => 'L', "CUC" => 'L', "CUG" => 'L', "CUU" => 'L',
            "GAA" => 'E', "GAC" => 'D', "GAG" => 'E', "GAU" => 'D',
            "GCA" => 'A', "GCC" => 'A', "GCG" => 'A', "GCU" => 'A',
            "GGA" => 'G', "GGC" => 'G', "GGG" => 'G', "GGU" => 'G',
            "GUA" => 'V', "GUC" => 'V', "GUG" => 'V', "GUU" => 'V',
            "UAA" => '*', "UAC" => 'Y', "UAG" => '*', "UAU" => 'Y',
            "UCA" => 'S', "UCC" => 'S', "UCG" => 'S', "UCU" => 'S',
            "UGA" => '*', "UGC" => 'C', "UGG" => 'W', "UGU" => 'C',
            "UUA" => 'L', "UUC" => 'F', "UUG" => 'L', "UUU" => 'F',
        ]

        function string_translate(seq::String)
            @assert length(seq) % 3 == 0
            aaseq = Array(Char, div(length(seq), 3))
            for i in 1:3:length(seq) - 3 + 1
                aaseq[div(i, 3) + 1] = standard_genetic_code_dict[seq[i:i+2]]
            end
            return convert(String, aaseq)
        end

        function check_translate(seq::String)
            return string_translate(seq) == convert(String, translate(RNASequence(seq)))
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_translate(random_translatable_rna(len)) for _ in 1:reps]) => true
        end

        @fact_throws translate(dna"ACGTACGTA") # can't translate DNA
        @fact_throws translate(rna"ACGUACGU") # can't translate non-multiples of three
        @fact_throws translate(rna"ACGUACGUN") # can't translate N
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

    context("Counting") do
        function string_nucleotide_count(::Type{DNANucleotide}, seq::String)
            counts = [
                DNA_A => 0,
                DNA_C => 0,
                DNA_G => 0,
                DNA_T => 0,
                DNA_N => 0 ]
            for c in seq
                counts[convert(DNANucleotide, c)] += 1
            end

            return counts
        end

        function string_nucleotide_count(::Type{RNANucleotide}, seq::String)
            counts = [
                RNA_A => 0,
                RNA_C => 0,
                RNA_G => 0,
                RNA_U => 0,
                RNA_N => 0 ]
            for c in seq
                counts[convert(RNANucleotide, c)] += 1
            end

            return counts
        end

        function check_nucleotide_count(::Type{DNANucleotide}, seq::String)
            string_counts = string_nucleotide_count(DNANucleotide, seq)
            seq_counts = nucleotide_count(DNASequence(seq))
            return string_counts[DNA_A] == seq_counts[DNA_A] &&
                   string_counts[DNA_C] == seq_counts[DNA_C] &&
                   string_counts[DNA_G] == seq_counts[DNA_G] &&
                   string_counts[DNA_T] == seq_counts[DNA_T] &&
                   string_counts[DNA_N] == seq_counts[DNA_N]
        end

        function check_nucleotide_count(::Type{RNANucleotide}, seq::String)
            string_counts = string_nucleotide_count(RNANucleotide, seq)
            seq_counts = nucleotide_count(RNASequence(seq))
            return string_counts[RNA_A] == seq_counts[RNA_A] &&
                   string_counts[RNA_C] == seq_counts[RNA_C] &&
                   string_counts[RNA_G] == seq_counts[RNA_G] &&
                   string_counts[RNA_U] == seq_counts[RNA_U] &&
                   string_counts[RNA_N] == seq_counts[RNA_N]
        end

        reps = 10
        for len in [1, 10, 32, 1000, 10000, 100000]
            @fact all([check_nucleotide_count(DNANucleotide, random_dna(len))
                       for _ in 1:reps]) => true
            @fact all([check_nucleotide_count(RNANucleotide, random_rna(len))
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
        function string_eachkmer(seq::String, k, step=1)
            kmers = String[]
            i = 1
            for i in 1:step:length(seq) - k + 1
                subseq = seq[i:i + k - 1]
                if !in('N', subseq)
                    push!(kmers, subseq)
                end
            end
            return kmers
        end

        function check_eachkmer(T, seq::String, k, step=1)
            xs = [convert(String, x) for x in collect(eachkmer(NucleotideSequence{T}(seq), k, step))]
            ys = string_eachkmer(seq, k, step)
            return xs == ys
        end

        reps = 10
        len = 10000
        for k in [1, 3, 16, 32]
            @fact all([check_eachkmer(DNANucleotide, random_dna(len), k) for _ in 1:reps]) => true
            @fact all([check_eachkmer(RNANucleotide, random_rna(len), k) for _ in 1:reps]) => true
        end

        len = 10000
        for k in [1, 3, 16, 32]
            @fact all([check_eachkmer(DNANucleotide, random_dna(len), k, 3) for _ in 1:reps]) => true
            @fact all([check_eachkmer(RNANucleotide, random_rna(len), k, 3) for _ in 1:reps]) => true
        end

        @fact isempty(collect(eachkmer(dna"", 1))) => true
        @fact isempty(collect(eachkmer(dna"NNNNNNNNNN", 1))) => true
        @fact isempty(collect(eachkmer(dna"ACGT", 0))) => true
        @fact_throws eachkmer(dna"ACGT", -1)
        @fact_throws eachkmer(dna"ACGT", 33)
    end
end



end # TestSeq
