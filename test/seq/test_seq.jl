module TestSeq

using Compat
using FactCheck
using YAML
using Bio
using Bio.Seq
import ..get_bio_fmt_specimens

# Return a random DNA/RNA sequence of the given length
function random_seq(n::Integer, nts, probs)
    cumprobs = cumsum(probs)
    x = Array(Char, n)
    for i in 1:n
        x[i] = nts[searchsorted(cumprobs, rand()).start]
    end
    return convert(String, x)
end


function random_array(n::Integer, elements, probs)
    cumprobs = cumsum(probs)
    x = Array(eltype(elements), n)
    for i in 1:n
        x[i] = elements[searchsorted(cumprobs, rand()).start]
    end
    return x
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
        "UUU",
        # translatable ambiguities in the standard code
        "CUN", "CCN", "CGN", "ACN",
        "GUN", "GCN", "GGN", "UCN"
        ]


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


function random_dna_kmer_nucleotides(len)
    return random_array(len, [DNA_A, DNA_C, DNA_G, DNA_T],
                        [0.25, 0.25, 0.25, 0.25])
end


function random_rna_kmer_nucleotides(len)
    return random_array(len, [RNA_A, RNA_C, RNA_G, RNA_U],
                        [0.25, 0.25, 0.25, 0.25])
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

facts("Nucleotides") do
    context("Conversions") do
        context("Uint8") do
            context("DNA conversions from Uint8") do
                @fact convert(DNANucleotide, @compat UInt8(0)) => DNA_A
                @fact convert(DNANucleotide, @compat UInt8(1)) => DNA_C
                @fact convert(DNANucleotide, @compat UInt8(2)) => DNA_G
                @fact convert(DNANucleotide, @compat UInt8(3)) => DNA_T
                @fact convert(DNANucleotide, @compat UInt8(4)) => DNA_N
            end

            context("RNA conversions from Uint8") do
                @fact convert(RNANucleotide, @compat UInt8(0)) => RNA_A
                @fact convert(RNANucleotide, @compat UInt8(1)) => RNA_C
                @fact convert(RNANucleotide, @compat UInt8(2)) => RNA_G
                @fact convert(RNANucleotide, @compat UInt8(3)) => RNA_U
                @fact convert(RNANucleotide, @compat UInt8(4)) => RNA_N
            end

            context("DNA conversions to Uint8") do
                @fact convert(Uint8, DNA_A) => @compat UInt8(0)
                @fact convert(Uint8, DNA_C) => @compat UInt8(1)
                @fact convert(Uint8, DNA_G) => @compat UInt8(2)
                @fact convert(Uint8, DNA_T) => @compat UInt8(3)
                @fact convert(Uint8, DNA_N) => @compat UInt8(4)
            end

            context("RNA conversions to Uint8") do
                @fact convert(Uint8, RNA_A) => @compat UInt8(0)
                @fact convert(Uint8, RNA_C) => @compat UInt8(1)
                @fact convert(Uint8, RNA_G) => @compat UInt8(2)
                @fact convert(Uint8, RNA_U) => @compat UInt8(3)
                @fact convert(Uint8, RNA_N) => @compat UInt8(4)
            end
        end

        context("Uint64") do
            context("DNA conversions from Uint64") do
                @fact convert(DNANucleotide, @compat UInt64(0)) => DNA_A
                @fact convert(DNANucleotide, @compat UInt64(1)) => DNA_C
                @fact convert(DNANucleotide, @compat UInt64(2)) => DNA_G
                @fact convert(DNANucleotide, @compat UInt64(3)) => DNA_T
                @fact convert(DNANucleotide, @compat UInt64(4)) => DNA_N
            end

            context("RNA conversions from Uint64") do
                @fact convert(RNANucleotide, @compat UInt64(0)) => RNA_A
                @fact convert(RNANucleotide, @compat UInt64(1)) => RNA_C
                @fact convert(RNANucleotide, @compat UInt64(2)) => RNA_G
                @fact convert(RNANucleotide, @compat UInt64(3)) => RNA_U
                @fact convert(RNANucleotide, @compat UInt64(4)) => RNA_N
            end

            context("DNA conversions to Uint64") do
                @fact convert(Uint64, DNA_A) => @compat UInt64(0)
                @fact convert(Uint64, DNA_C) => @compat UInt64(1)
                @fact convert(Uint64, DNA_G) => @compat UInt64(2)
                @fact convert(Uint64, DNA_T) => @compat UInt64(3)
                @fact convert(Uint64, DNA_N) => @compat UInt64(4)
            end

            context("RNA conversions to Uint64") do
                @fact convert(Uint64, RNA_A) => @compat UInt64(0)
                @fact convert(Uint64, RNA_C) => @compat UInt64(1)
                @fact convert(Uint64, RNA_G) => @compat UInt64(2)
                @fact convert(Uint64, RNA_U) => @compat UInt64(3)
                @fact convert(Uint64, RNA_N) => @compat UInt64(4)
            end
        end

        context("Char") do
            context("DNA conversions from Char") do
                @fact convert(DNANucleotide, 'A') => DNA_A
                @fact convert(DNANucleotide, 'C') => DNA_C
                @fact convert(DNANucleotide, 'G') => DNA_G
                @fact convert(DNANucleotide, 'T') => DNA_T
                @fact convert(DNANucleotide, 'N') => DNA_N
            end

            context("RNA conversions from Char") do
                @fact convert(RNANucleotide, 'A') => RNA_A
                @fact convert(RNANucleotide, 'C') => RNA_C
                @fact convert(RNANucleotide, 'G') => RNA_G
                @fact convert(RNANucleotide, 'U') => RNA_U
                @fact convert(RNANucleotide, 'N') => RNA_N
            end

            context("DNA conversions to Char") do
                @fact convert(Char, DNA_A) => 'A'
                @fact convert(Char, DNA_C) => 'C'
                @fact convert(Char, DNA_G) => 'G'
                @fact convert(Char, DNA_T) => 'T'
                @fact convert(Char, DNA_N) => 'N'
            end

            context("DNA conversions to Char") do
                @fact convert(Char, RNA_A) => 'A'
                @fact convert(Char, RNA_C) => 'C'
                @fact convert(Char, RNA_G) => 'G'
                @fact convert(Char, RNA_U) => 'U'
                @fact convert(Char, RNA_N) => 'N'
            end
        end
    end

    context("Sequences") do
        function dna_complement(seq::String)
            seqc = Array(Char, length(seq))
            for (i, c) in enumerate(seq)
                if c     ==   'A'
                    seqc[i] = 'T'
                elseif c ==   'C'
                    seqc[i] = 'G'
                elseif c ==   'G'
                    seqc[i] = 'C'
                elseif c ==   'T'
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

        function check_reversal(T, seq)
            return reverse(seq) == convert(String, reverse(convert(T, seq)))
        end

        function check_dna_complement(T, seq)
            return dna_complement(seq) ==
                convert(String, complement(convert(T, seq)))
        end

        function check_rna_complement(T, seq)
            return rna_complement(seq) ==
                convert(String, complement(convert(T, seq)))
        end

        function check_dna_revcomp(T, seq)
            return reverse(dna_complement(seq)) ==
                convert(String, reverse_complement(convert(T, seq)))
        end

        function check_rna_revcomp(T, seq)
            return reverse(rna_complement(seq)) ==
                convert(String, reverse_complement(convert(T, seq)))
        end

        function check_mismatches(T, a, b)
            count = 0
            for (ca, cb) in zip(a, b)
                if ca != cb
                    count += 1
                end
            end
            return mismatches(convert(T, a), convert(T, b)) == count
        end

        context("Nucleotide Sequences") do
            reps = 10
            context("Construction and Conversions") do
                context("Constructing empty sequences") do
                    # Check construction of empty nucleotide sequences
                    #  using RNASequence and DNASequence functions
                    @fact RNASequence() => NucleotideSequence(RNANucleotide)
                    @fact DNASequence() => NucleotideSequence(DNANucleotide)
                end

                context("Conversion from/to Strings") do
                    # Check that sequences in strings survive round trip conversion:
                    #   String → NucleotideSequence → String
                    function check_string_construction(T::Type, seq::String)
                        return convert(String, NucleotideSequence{T}(seq)) == uppercase(seq)
                    end

                    for len in [0, 1, 10, 32, 1000, 10000, 100000]
                        @fact all([check_string_construction(DNANucleotide, random_dna(len)) for _ in 1:reps]) => true
                        @fact all([check_string_construction(RNANucleotide, random_rna(len)) for _ in 1:reps]) => true
                        @fact all([check_string_construction(DNANucleotide, lowercase(random_dna(len))) for _ in 1:reps]) => true
                        @fact all([check_string_construction(RNANucleotide, lowercase(random_rna(len))) for _ in 1:reps]) => true
                    end

                    # Non-nucleotide characters should throw
                    @fact_throws DNASequence("ACCNNCATTTTTTAGATXATAG")
                    @fact_throws RNASequence("ACCNNCATTTTTTAGATXATAG")
                end

                context("Conversion between RNA and DNA") do
                    @fact convert(RNASequence, DNASequence("ACGTN")) => rna"ACGUN"
                    @fact convert(DNASequence, RNASequence("ACGUN")) => dna"ACGTN"
                end

                context("Concatenation") do
                    function check_concatenation(::Type{DNANucleotide}, n)
                        chunks = [random_dna(rand(100:300)) for i in 1:n]
                        parts = Any[]
                        for i in 1:n
                            start = rand(1:length(chunks[i]))
                            stop = rand(start:length(chunks[i]))
                            push!(parts, start:stop)
                        end

                        str = string([chunk[parts[i]]
                                      for (i, chunk) in enumerate(chunks)]...)

                        seq = *([DNASequence(chunk)[parts[i]]
                                 for (i, chunk) in enumerate(chunks)]...)

                        return convert(String, seq) == uppercase(str)
                    end

                    @fact all([check_concatenation(DNANucleotide, rand(1:10)) for _ in 1:100]) => true
                end

                context("Repetition") do
                    function check_repetition(::Type{DNANucleotide}, n)
                        chunk = random_dna(rand(100:300))
                        start = rand(1:length(chunk))
                        stop = rand(start:length(chunk))

                        str = chunk[start:stop] ^ n
                        seq = DNASequence(chunk)[start:stop] ^ n
                        return convert(String, seq) == uppercase(str)
                    end

                    @fact all([check_repetition(DNANucleotide, rand(1:10)) for _ in 1:100]) => true
                end
            end

            context("Equality") do
                reps = 10
                function check_seq_equality(len)
                    a = random_dna(len)
                    return a == copy(a)
                end

                for len in [1, 10, 32, 1000]
                    @fact all([check_seq_equality(len) for _ in 1:reps]) => true
                end

                a = b = dna"ACTGN"
                @fact ==(a, b) => true
                @fact ==(dna"ACTGN", dna"ACTGN") => true
                @fact ==(dna"ACTGN", dna"ACTGA") => false
                @fact ==(dna"ACTGN", dna"ACTG") => false
                @fact ==(dna"ACTG", dna"ACTGN") => false

                c = d = rna"ACUGN"
                @fact ==(c, d) => true
                @fact ==(rna"ACUGN", rna"ACUGN") => true
                @fact ==(rna"ACUGN", rna"ACUGA") => false
                @fact ==(rna"ACUGN", rna"ACUG") => false
                @fact ==(rna"ACUG", rna"ACUGN") => false
            end

            context("Length") do
                for len in [0, 1, 10, 32, 1000, 10000, 100000]
                    @fact length(DNASequence(random_dna(len))) => len
                    @fact endof(DNASequence(random_dna(len))) => len
                    @fact length(DNASequence(random_dna(len))) => len
                    @fact endof(DNASequence(random_dna(len))) => len
                end
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

            context("Access and Iterations") do
                dna_seq = dna"ACTG"
                rna_seq = rna"ACUG"

                context("Access DNA Sequence") do
                    # Access indexes out of bounds
                    @fact_throws dna_seq[-1]
                    @fact_throws dna_seq[0]
                    @fact_throws dna_seq[5]
                    @fact_throws getindex(dna_seq,-1)
                    @fact_throws getindex(dna_seq, 0)
                    @fact_throws getindex(dna_seq, 5)
                end

                context("Iteration through DNA Sequence") do
                    @fact start(dna"ACNTG") => (1,3)
                    @fact start(dna"")      => (1,1)

                    @fact next(dna"ACTGN", (2,5)) => (DNA_C, (3,5))
                    @fact next(dna"ACTGN", (5,5)) => (DNA_N, (6,6))

                    @fact done(dna"", (1,1))       => true
                    @fact done(dna"ACTGN", (2,5))  => false
                    @fact done(dna"ACTGN", (5,5))  => false
                    @fact done(dna"ACTGN", (6,5))  => true
                    @fact done(dna"ACTGN", (0,5))  => false

                    dna_vector = [DNA_A, DNA_C, DNA_T, DNA_G]
                    @fact all([nucleotide == dna_vector[i] for (i, nucleotide) in enumerate(dna_seq)]) =>  true
                end

                context("Access RNA Sequence") do
                    # Access indexes out of bounds
                    @fact_throws rna_seq[-1]
                    @fact_throws rna_seq[0]
                    @fact_throws rna_seq[5]
                    @fact_throws getindex(rna_seq, -1)
                    @fact_throws getindex(rna_seq, 0)
                    @fact_throws getindex(rna_seq, 5)
                end

                context("Iteration through RNA Sequence") do
                    @fact start(rna"ACNUG") => (1,3)
                    @fact start(rna"")      => (1,1)

                    @fact next(rna"ACUGN", (2,5)) => (RNA_C, (3,5))
                    @fact next(rna"ACUGN", (5,5)) => (RNA_N, (6,6))

                    @fact done(rna"", (1,1))       => true
                    @fact done(rna"ACUGN", (2,5))  => false
                    @fact done(rna"ACUGN", (5,5))  => false
                    @fact done(rna"ACUGN", (6,5))  => true
                    @fact done(rna"ACUGN", (0,5))  => false

                    # Iteration through RNA Sequence

                    rna_vector = [RNA_A, RNA_C, RNA_U, RNA_G]
                    @fact all([nucleotide == rna_vector[i] for (i, nucleotide) in enumerate(rna_seq)]) =>  true
                end

                context("Indexing with Ranges") do
                    @fact getindex(dna"ACTGNACTGN", 1:5) => dna"ACTGN"
                    @fact getindex(rna"ACUGNACUGN", 1:5) => rna"ACUGN"
                    @fact getindex(dna"ACTGNACTGN", 5:1) => dna""
                    @fact getindex(rna"ACUGNACUGN", 5:1) => rna""
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

                context("Subsequence Construction from Ranges") do
                    # Subsequence from range
                    @fact RNASequence(rna"AUCGAUCG", 5:8) == RNASequence("AUCG") => true
                    @fact DNASequence(dna"ATCGATCG", 5:8) == DNASequence("ATCG") => true

                    # Invalid ranges
                    @fact_throws RNASequence(rna"AUCGAUCG", 5:10)
                    @fact_throws DNASequence(dna"ATCGATCG", 5:10)

                    # Empty ranges
                    @fact RNASequence(rna"AUCGAUCG", 5:4) == RNASequence() => true
                    @fact DNASequence(dna"ATCGATCG", 5:4) == DNASequence() => true
                end
            end

            context("Transformations") do
                context("Reversal") do
                    for len in [0, 1, 10, 32, 1000, 10000, 100000]
                        @fact all([check_reversal(DNASequence, random_dna(len)) for _ in 1:reps]) => true
                        @fact all([check_reversal(RNASequence, random_rna(len)) for _ in 1:reps]) => true
                    end
                end

                context("Complement") do
                    for len in [1, 10, 32, 1000, 10000, 100000]
                        @fact all([check_dna_complement(DNASequence, random_dna(len)) for _ in 1:reps]) => true
                        @fact all([check_rna_complement(RNASequence, random_rna(len)) for _ in 1:reps]) => true
                    end
                end

                context("Reverse Complement") do
                    for len in [1, 10, 32, 1000, 10000, 100000]
                        @fact all([check_dna_revcomp(DNASequence, random_dna(len)) for _ in 1:reps]) => true
                        @fact all([check_rna_revcomp(RNASequence, random_rna(len)) for _ in 1:reps]) => true
                    end
                end
            end

            context("Mismatches") do
                for len in [1, 10, 32, 1000, 10000, 100000]
                    @fact all([check_mismatches(DNASequence, random_dna(len), random_dna(len))
                               for _ in 1:reps]) => true
                    @fact all([check_mismatches(RNASequence, random_rna(len), random_rna(len))
                               for _ in 1:reps]) => true
                end
            end
        end

        context("SequenceNIterator") do
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

            dna_seq   = dna"ANANANA"
            dna_niter = npositions(dna_seq)
            @fact start(dna_niter) => 2
            @fact next(dna_niter, 2) => (2,4)
            @fact done(dna_niter, 6) => false
            @fact done(dna_niter, 8) => true


            ns = [2,4,6]
            for (i, n) in enumerate(dna_niter)
                @fact n => ns[i]
            end

            rna_seq   = rna"ANANANA"
            rna_niter = npositions(rna_seq)
            @fact start(rna_niter) => 2
            @fact next(rna_niter, 2) => (2,4)
            @fact done(rna_niter, 6) => false
            @fact done(rna_niter, 8) => true

            ns = [2,4,6]
            for (i, n) in enumerate(rna_niter)
                @fact n => ns[i]
            end

        end

        context("Kmer") do
            reps = 100
            context("Construction and Conversions") do
                # Check that kmers in strings survive round trip conversion:
                #   Uint64 → Kmer → Uint64
                function check_uint64_convertion(T::Type, n::Uint64, len::Int)
                    return convert(Uint64, convert(Kmer{T, len}, n)) === n
                end

                # Check that kmers in strings survive round trip conversion:
                #   String → Kmer → String
                function check_string_construction(T::Type, seq::String)
                    return convert(String, convert(Kmer{T}, seq)) == uppercase(seq)
                end

                # Check that dnakmers can be constructed from a DNASequence
                #   DNASequence → Kmer → DNASequence
                function check_dnasequence_construction(seq::DNASequence)
                    return convert(DNASequence, convert(DNAKmer, seq)) == seq
                end

                # Check that rnakmers can be constructed from a RNASequence
                #   RNASequence → Kmer → RNASequence
                function check_rnasequence_construction(seq::RNASequence)
                    return convert(RNASequence, convert(RNAKmer, seq)) == seq
                end

                # Check that kmers can be constructed from a NucleotideSequence
                #   NucleotideSequence → Kmer → NucleotideSequence
                function check_nucsequence_construction(seq::NucleotideSequence)
                    return convert(NucleotideSequence, convert(Kmer, seq)) == seq
                end

                # Check that kmers can be constructed from an array of nucleotides
                #   Vector{Nucleotide} → Kmer → Vector{Nucleotide}
                function check_nucarray_kmer{T <: Nucleotide}(seq::Vector{T})
                    return convert(String, [convert(Char, c) for c in seq]) ==
                           convert(String, kmer(seq...))
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
                    # Uint64 conversions
                    @fact all([check_uint64_convertion(DNANucleotide, rand(Uint64), len)
                              for _ in 1:reps]) => true
                    @fact all([check_uint64_convertion(RNANucleotide, rand(Uint64), len)
                              for _ in 1:reps]) => true

                    # String construction
                    @fact all([check_string_construction(DNANucleotide, random_dna_kmer(len))
                               for _ in 1:reps]) => true
                    @fact all([check_string_construction(RNANucleotide, random_rna_kmer(len))
                               for _ in 1:reps]) => true

                    # DNA/RNASequence Constructions
                    @fact all([check_dnasequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps]) => true
                    @fact all([check_rnasequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps]) => true

                    # NucleotideSequence Construction
                    @fact all([check_nucsequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps]) => true
                    @fact all([check_nucsequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps]) => true

                    # Construction from nucleotide arrays
                    if len > 0
                        @fact all([check_nucarray_kmer(random_dna_kmer_nucleotides(len))
                                   for _ in 1:reps]) => true
                        @fact all([check_nucarray_kmer(random_rna_kmer_nucleotides(len))
                                   for _ in 1:reps]) => true
                    end

                    # Roundabout conversions
                    @fact all([check_roundabout_construction(DNANucleotide, random_dna_kmer(len))
                               for _ in 1:reps]) => true
                    @fact all([check_roundabout_construction(RNANucleotide, random_rna_kmer(len))
                               for _ in 1:reps]) => true
                end

                @fact_throws kmer() # can't construct 0-mer using `kmer()`
                @fact_throws kmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
                @fact_throws kmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
                @fact_throws kmer(rna"ACGNU")# no Ns in kmers
                @fact_throws rnmakmer(rna"ACGNU")# no Ns in kmers
                @fact_throws kmer(dna"ACGNT") # no Ns in kmers
                @fact_throws dnmakmer(dna"ACGNT") # no Ns in kmers
                @fact_throws kmer(RNA_A, DNA_A) # no mixing of RNA and DNA
                @fact_throws kmer(random_rna(33)) # no kmer larger than 32nt
                @fact_throws kmer(random_dna(33)) # no kmer larger than 32nt
                @fact_throws kmer(RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 32nt
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U)
                @fact_throws kmer(DNA_A, DNA_C, DNA_G, DNA_T, # no kmer larger than 32nt
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T)

                context("From strings") do
                    @fact dnakmer("ACTG") => convert(Kmer, DNASequence("ACTG"))
                    @fact rnakmer("ACUG") => convert(Kmer, RNASequence("ACUG"))

                    # N is not allowed in Kmers
                    @fact_throws dnakmer("ACGTNACGT")
                    @fact_throws rnakmer("ACGUNACGU")
                end
            end

            context("Comparisons") do
                context("Equality") do
                    function check_seq_kmer_equality(len)
                        a = dnakmer(random_dna_kmer(len))
                        b = convert(DNASequence, a)
                        return a == b && b == a
                    end

                    for len in [1, 10, 32]
                        @fact all([check_seq_kmer_equality(len) for _ in 1:reps]) => true
                    end

                    # True negatives
                    @fact dnakmer("ACG") == rnakmer("ACG") => false
                    @fact dnakmer("T")   == rnakmer("U")   => false
                    @fact dnakmer("AC")  == dnakmer("AG")  => false
                    @fact rnakmer("AC")  == rnakmer("AG")  => false

                    @fact dnakmer("ACG") == rna"ACG" => false
                    @fact dnakmer("T")   == rna"U"   => false
                    @fact dnakmer("AC")  == dna"AG"  => false
                    @fact rnakmer("AC")  == rna"AG"  => false

                    @fact rna"ACG" == dnakmer("ACG") => false
                    @fact rna"U"   == dnakmer("T")   => false
                    @fact dna"AG"  == dnakmer("AC")  => false
                    @fact rna"AG"  == rnakmer("AC")  => false
                end

                context("Inequality") do
                    for len in [1, 10, 32]
                        @fact isless(convert(DNAKmer{1}, @compat UInt64(0)), convert(DNAKmer{1}, @compat UInt64(1))) => true
                        @fact isless(convert(DNAKmer{1}, @compat UInt64(0)), convert(DNAKmer{1}, @compat UInt64(0))) => false
                        @fact isless(convert(DNAKmer{1}, @compat UInt64(1)), convert(DNAKmer{1}, @compat UInt64(0))) => false

                        @fact isless(convert(RNAKmer{1}, @compat UInt64(0)), convert(RNAKmer{1}, @compat UInt64(1))) => true
                        @fact isless(convert(RNAKmer{1}, @compat UInt64(0)), convert(RNAKmer{1}, @compat UInt64(0))) => false
                        @fact isless(convert(RNAKmer{1}, @compat UInt64(1)), convert(RNAKmer{1}, @compat UInt64(0))) => false
                    end
                end
            end

            context("Length") do
                for len in [0, 1, 16, 32]
                    @fact length(dnakmer(random_dna_kmer(len))) => len
                    @fact length(rnakmer(random_rna_kmer(len))) => len
                end
            end

            context("Access and Iterations") do
                dna_kmer = dnakmer("ACTG")
                rna_kmer = rnakmer("ACUG")

                context("Access DNA Kmer") do
                    @fact dna_kmer[1] => DNA_A
                    @fact dna_kmer[2] => DNA_C
                    @fact dna_kmer[3] => DNA_T
                    @fact dna_kmer[4] => DNA_G

                    # Access indexes out of bounds
                    @fact_throws dna_kmer[-1]
                    @fact_throws dna_kmer[0]
                    @fact_throws dna_kmer[5]
                    @fact_throws getindex(dna_kmer,-1)
                    @fact_throws getindex(dna_kmer, 0)
                    @fact_throws getindex(dna_kmer, 5)
                end

                context("Iteration through DNA Kmer") do
                    @fact start(dnakmer("ACTG"))  => 1
                    @fact start(dnakmer(""))      => 1

                    @fact next(dnakmer("ACTG"), 1) => (DNA_A, 2)
                    @fact next(dnakmer("ACTG"), 4) => (DNA_G, 5)

                    @fact done(dnakmer(""), 1)      => true
                    @fact done(dnakmer("ACTG"), 1)  => false
                    @fact done(dnakmer("ACTG"), 4)  => false
                    @fact done(dnakmer("ACTG"), 5)  => true
                    @fact done(dnakmer("ACTG"), -1) => false


                    dna_kmer_vector = [DNA_A, DNA_C, DNA_T, DNA_G]
                    @fact all([nucleotide == dna_kmer_vector[i] for (i, nucleotide) in enumerate(dna_kmer)]) =>  true
                end

                context("Access RNA Kmer") do
                    @fact rna_kmer[1] => RNA_A
                    @fact rna_kmer[2] => RNA_C
                    @fact rna_kmer[3] => RNA_U
                    @fact rna_kmer[4] => RNA_G

                    # Access indexes out of bounds
                    @fact_throws rna_kmer[-1]
                    @fact_throws rna_kmer[0]
                    @fact_throws rna_kmer[5]
                    @fact_throws getindex(rna_kmer, -1)
                    @fact_throws getindex(rna_kmer, 0)
                    @fact_throws getindex(rna_kmer, 5)
                end

                context("Iteration through RNA Kmer") do
                    @fact start(rnakmer("ACUG"))  => 1
                    @fact start(rnakmer(""))      => 1

                    @fact next(rnakmer("ACUG"), 1) => (RNA_A, 2)
                    @fact next(rnakmer("ACUG"), 4) => (RNA_G, 5)

                    @fact done(rnakmer(""), 1)      => true
                    @fact done(rnakmer("ACUG"), 1)  => false
                    @fact done(rnakmer("ACUG"), 4)  => false
                    @fact done(rnakmer("ACUG"), 5)  => true
                    @fact done(rnakmer("ACUG"), -1) => false

                    rna_kmer_vector = [RNA_A, RNA_C, RNA_U, RNA_G]
                    @fact all([nucleotide == rna_kmer_vector[i] for (i, nucleotide) in enumerate(rna_kmer)]) =>  true
                end
            end

            reps = 10
            context("Transformations") do
                context("Reversal") do
                    for len in [0, 1, 16, 32]
                        @fact all([check_reversal(DNAKmer, random_dna_kmer(len)) for _ in 1:reps]) => true
                        @fact all([check_reversal(RNAKmer, random_rna_kmer(len)) for _ in 1:reps]) => true
                    end
                end

                context("Complement") do
                    for len in [0, 1, 16, 32]
                        @fact all([check_dna_complement(DNAKmer, random_dna_kmer(len)) for _ in 1:reps]) => true
                        @fact all([check_rna_complement(RNAKmer, random_rna_kmer(len)) for _ in 1:reps]) => true
                    end
                end

                context("Reverse Complement") do
                    for len in [0, 1, 16, 32]
                        @fact all([check_dna_revcomp(DNAKmer, random_dna_kmer(len)) for _ in 1:reps]) => true
                        @fact all([check_rna_revcomp(RNAKmer, random_rna_kmer(len)) for _ in 1:reps]) => true
                    end
                end
            end

            context("Mismatches") do
                for len in [0, 1, 16, 32]
                    @fact all([check_mismatches(DNAKmer, random_dna_kmer(len), random_dna_kmer(len))
                              for _ in 1:reps]) => true
                    @fact all([check_mismatches(RNAKmer, random_rna_kmer(len), random_rna_kmer(len))
                              for _ in 1:reps]) => true
                end
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
                xs = [convert(String, x) for (i, x) in collect(eachkmer(NucleotideSequence{T}(seq), k, step))]
                ys = string_eachkmer(seq, k, step)
                return xs == ys
            end

            reps = 10
            len = 10000

            for k in [1, 3, 16, 32]
                @fact all([check_eachkmer(DNANucleotide, random_dna(len), k) for _ in 1:reps]) => true
                @fact all([check_eachkmer(RNANucleotide, random_rna(len), k) for _ in 1:reps]) => true
            end

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

        context("Counting") do
            function string_nucleotide_count(::Type{DNANucleotide}, seq::String)
                counts = @compat Dict{DNANucleotide, Int}(
                    DNA_A => 0,
                    DNA_C => 0,
                    DNA_G => 0,
                    DNA_T => 0,
                    DNA_N => 0 )
                for c in seq
                    counts[convert(DNANucleotide, c)] += 1
                end

                return counts
            end

            function string_nucleotide_count(::Type{RNANucleotide}, seq::String)
                counts = @compat Dict{RNANucleotide, Int}(
                    RNA_A => 0,
                    RNA_C => 0,
                    RNA_G => 0,
                    RNA_U => 0,
                    RNA_N => 0 )
                for c in seq
                    counts[convert(RNANucleotide, c)] += 1
                end

                return counts
            end

            function check_nucleotide_count(::Type{DNANucleotide}, seq::String)
                string_counts = string_nucleotide_count(DNANucleotide, seq)
                seq_counts = NucleotideCounts(DNASequence(seq))
                return string_counts[DNA_A] == seq_counts[DNA_A] &&
                       string_counts[DNA_C] == seq_counts[DNA_C] &&
                       string_counts[DNA_G] == seq_counts[DNA_G] &&
                       string_counts[DNA_T] == seq_counts[DNA_T] &&
                       string_counts[DNA_N] == seq_counts[DNA_N]
            end

            function check_nucleotide_count(::Type{RNANucleotide}, seq::String)
                string_counts = string_nucleotide_count(RNANucleotide, seq)
                seq_counts = NucleotideCounts(RNASequence(seq))
                return string_counts[RNA_A] == seq_counts[RNA_A] &&
                       string_counts[RNA_C] == seq_counts[RNA_C] &&
                       string_counts[RNA_G] == seq_counts[RNA_G] &&
                       string_counts[RNA_U] == seq_counts[RNA_U] &&
                       string_counts[RNA_N] == seq_counts[RNA_N]
            end

            function check_kmer_nucleotide_count(::Type{DNANucleotide}, seq::String)
                string_counts = string_nucleotide_count(DNANucleotide, seq)
                kmer_counts = NucleotideCounts(dnakmer(seq))
                return string_counts[DNA_A] == kmer_counts[DNA_A] &&
                       string_counts[DNA_C] == kmer_counts[DNA_C] &&
                       string_counts[DNA_G] == kmer_counts[DNA_G] &&
                       string_counts[DNA_T] == kmer_counts[DNA_T] &&
                       string_counts[DNA_N] == kmer_counts[DNA_N]
            end

            function check_kmer_nucleotide_count(::Type{RNANucleotide}, seq::String)
                string_counts = string_nucleotide_count(RNANucleotide, seq)
                kmer_counts = NucleotideCounts(rnakmer(seq))
                return string_counts[RNA_A] == kmer_counts[RNA_A] &&
                       string_counts[RNA_C] == kmer_counts[RNA_C] &&
                       string_counts[RNA_G] == kmer_counts[RNA_G] &&
                       string_counts[RNA_U] == kmer_counts[RNA_U] &&
                       string_counts[RNA_N] == kmer_counts[RNA_N]
            end

            reps = 10
            for len in [1, 10, 32, 1000, 10000, 100000]
                @fact all([check_nucleotide_count(DNANucleotide, random_dna(len))
                           for _ in 1:reps]) => true
                @fact all([check_nucleotide_count(RNANucleotide, random_rna(len))
                           for _ in 1:reps]) => true
            end

            for len in [1, 10, 32]
                @fact all([check_kmer_nucleotide_count(DNANucleotide, random_dna_kmer(len))
                           for _ in 1:reps]) => true
                @fact all([check_kmer_nucleotide_count(RNANucleotide, random_rna_kmer(len))
                           for _ in 1:reps]) => true
            end
        end
    end
end

facts("Aminoacids") do
    context("AminoAcid Sequences") do
        reps = 10
        context("Construction") do
            # Non-aa characters should throw
            @fact_throws AminoAcidSequence("ATGHLMYZZACAGNM")

            # Check that sequences in strings survive round trip conversion:
            #   String → AminoAcidSequence → String
            function check_string_construction(seq::String)
                return convert(String, AminoAcidSequence(seq)) == uppercase(seq)
            end

            for len in [0, 1, 10, 32, 1000, 10000, 100000]
                @fact all([check_string_construction(random_aa(len)) for _ in 1:reps]) => true
                @fact all([check_string_construction(lowercase(random_aa(len))) for _ in 1:reps]) => true
            end

            # Check creation of empty
        end

        context("Concatenation") do
            function check_concatenation(n)
                chunks = [random_aa(rand(100:300)) for i in 1:n]
                parts = Any[]
                for i in 1:n
                    start = rand(1:length(chunks[i]))
                    stop = rand(start:length(chunks[i]))
                    push!(parts, start:stop)
                end

                str = string([chunk[parts[i]]
                              for (i, chunk) in enumerate(chunks)]...)

                seq = *([AminoAcidSequence(chunk)[parts[i]]
                         for (i, chunk) in enumerate(chunks)]...)

                return convert(String, seq) == uppercase(str)
            end

            @fact all([check_concatenation(rand(1:10)) for _ in 1:100]) => true
        end

        context("Repetition") do
            function check_repetition(n)
                chunk = random_aa(rand(100:300))
                start = rand(1:length(chunk))
                stop = rand(start:length(chunk))

                str = chunk[start:stop] ^ n
                seq = AminoAcidSequence(chunk)[start:stop] ^ n

                return convert(String, seq) == uppercase(str)
            end

            @fact all([check_repetition(rand(1:10)) for _ in 1:100]) => true
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

    context("Parsers") do
        context("Valid Cases") do
            # case-insensitive and ignores spaces
            @fact parse(AminoAcid, "a") => AA_A
            @fact parse(AminoAcid, "Ala") => AA_A
            @fact parse(AminoAcid, "aLa ") => AA_A
            @fact parse(AminoAcid, " alA ") => AA_A
            @fact parse(AminoAcid, "\tAlA\n") => AA_A
            @fact parse(AminoAcid, "x") => AA_X
            @fact parse(AminoAcid, "X") => AA_X
            aas = [
                ("A", "ALA", AA_A),
                ("R", "ARG", AA_R),
                ("N", "ASN", AA_N),
                ("D", "ASP", AA_D),
                ("C", "CYS", AA_C),
                ("E", "GLU", AA_E),
                ("Q", "GLN", AA_Q),
                ("G", "GLY", AA_G),
                ("H", "HIS", AA_H),
                ("I", "ILE", AA_I),
                ("L", "LEU", AA_L),
                ("K", "LYS", AA_K),
                ("M", "MET", AA_M),
                ("F", "PHE", AA_F),
                ("P", "PRO", AA_P),
                ("S", "SER", AA_S),
                ("T", "THR", AA_T),
                ("W", "TRP", AA_W),
                ("Y", "TYR", AA_Y),
                ("V", "VAL", AA_V),
            ]
            @fact length(aas) => 20
            for (one, three, aa) in aas
                @fact parse(AminoAcid, one) => aa
                @fact parse(AminoAcid, three) => aa
            end
        end

        context("Invalid Cases") do
            @fact_throws parse(AminoAcid, "")
            @fact_throws parse(AminoAcid, "AL")
            @fact_throws parse(AminoAcid, "LA")
            @fact_throws parse(AminoAcid, "ALAA")
        end
    end
end

facts("Translation") do
    # crummy string translation to test against
    standard_genetic_code_dict = @compat Dict{String, Char}(
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

        # translatable ambiguities in the standard code
        "CUN" => 'L', "CCN" => 'P', "CGN" => 'R', "ACN" => 'T',
        "GUN" => 'V', "GCN" => 'A', "GGN" => 'G', "UCN" => 'S'

    )

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
    #for len in [1, 10, 32, 1000, 10000, 100000]
    for len in [1, 10, 32]
        @fact all([check_translate(random_translatable_rna(len)) for _ in 1:reps]) => true
    end

    @fact_throws translate(dna"ACGTACGTA") # can't translate DNA
    @fact_throws translate(rna"ACGUACGU")  # can't translate non-multiples of three
    @fact_throws translate(rna"ACGUACGNU") # can't translate N
end


facts("Sequence Parsing") do
    context("FASTA Parsing") do
        get_bio_fmt_specimens()

        function check_fasta_parse(filename)
            # Reading from a stream
            for seqrec in read(open(filename), FASTA)
            end

            # Reading from a memory mapped file
            for seqrec in read(filename, FASTA, memory_map=true)
            end

            return true
        end

        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "FASTA")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            tags = specimen["tags"]
            # currently unsupported features
            if contains(tags, "gaps") || contains(tags, "comments") || contains(tags, "ambiguity")
                continue
            end
            @fact check_fasta_parse(joinpath(path, specimen["filename"])) => true
        end
    end

    context("FASTQ Parsing") do
        get_bio_fmt_specimens()

        function check_fastq_parse(filename)
            # Reading from a stream
            for seqrec in read(open(filename), FASTQ)
            end

            # Reading from a memory mapped file
            for seqrec in read(filename, FASTQ, memory_map=true)
            end

            return true
        end

        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "FASTQ")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            tags = get(specimen, "tags", "")
            valid = get(specimen, "valid", true)
            # currently unsupported features
            if contains(tags, "rna") || contains(tags, "gaps") ||
               contains(tags, "comments") || contains(tags, "ambiguity")
                continue
            end
            if valid
                @fact check_fastq_parse(joinpath(path, specimen["filename"])) => true
            else
                @fact_throws check_fastq_parse(joinpath(path, specimen["filename"]))
            end
        end
    end
end

end # TestSeq
