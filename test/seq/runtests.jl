module TestSeq

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Bio.Seq,
    YAML,
    TestFunctions


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
    x = Array(AbstractString, n)
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


function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
end

@testset "Nucleotides" begin
    @testset "Conversions" begin
        @testset "UInt8" begin
            @testset "DNA conversions from UInt8" begin
                @test convert(DNANucleotide, UInt8(0)) == DNA_A
                @test convert(DNANucleotide, UInt8(1)) == DNA_C
                @test convert(DNANucleotide, UInt8(2)) == DNA_G
                @test convert(DNANucleotide, UInt8(3)) == DNA_T
                @test convert(DNANucleotide, UInt8(4)) == DNA_N
            end

            @testset "RNA conversions from UInt8" begin
                @test convert(RNANucleotide, UInt8(0)) == RNA_A
                @test convert(RNANucleotide, UInt8(1)) == RNA_C
                @test convert(RNANucleotide, UInt8(2)) == RNA_G
                @test convert(RNANucleotide, UInt8(3)) == RNA_U
                @test convert(RNANucleotide, UInt8(4)) == RNA_N
            end

            @testset "DNA conversions to UInt8" begin
                @test convert(UInt8, DNA_A) == UInt8(0)
                @test convert(UInt8, DNA_C) == UInt8(1)
                @test convert(UInt8, DNA_G) == UInt8(2)
                @test convert(UInt8, DNA_T) == UInt8(3)
                @test convert(UInt8, DNA_N) == UInt8(4)
            end

            @testset "RNA conversions to UInt8" begin
                @test convert(UInt8, RNA_A) == UInt8(0)
                @test convert(UInt8, RNA_C) == UInt8(1)
                @test convert(UInt8, RNA_G) == UInt8(2)
                @test convert(UInt8, RNA_U) == UInt8(3)
                @test convert(UInt8, RNA_N) == UInt8(4)
            end
        end

        @testset "UInt64" begin
            @testset "DNA conversions from UInt64" begin
                @test convert(DNANucleotide, UInt64(0)) == DNA_A
                @test convert(DNANucleotide, UInt64(1)) == DNA_C
                @test convert(DNANucleotide, UInt64(2)) == DNA_G
                @test convert(DNANucleotide, UInt64(3)) == DNA_T
                @test convert(DNANucleotide, UInt64(4)) == DNA_N
            end

            @testset "RNA conversions from UInt64" begin
                @test convert(RNANucleotide, UInt64(0)) == RNA_A
                @test convert(RNANucleotide, UInt64(1)) == RNA_C
                @test convert(RNANucleotide, UInt64(2)) == RNA_G
                @test convert(RNANucleotide, UInt64(3)) == RNA_U
                @test convert(RNANucleotide, UInt64(4)) == RNA_N
            end

            @testset "DNA conversions to UInt64" begin
                @test convert(UInt64, DNA_A) == UInt64(0)
                @test convert(UInt64, DNA_C) == UInt64(1)
                @test convert(UInt64, DNA_G) == UInt64(2)
                @test convert(UInt64, DNA_T) == UInt64(3)
                @test convert(UInt64, DNA_N) == UInt64(4)
            end

            @testset "RNA conversions to UInt64" begin
                @test convert(UInt64, RNA_A) == UInt64(0)
                @test convert(UInt64, RNA_C) == UInt64(1)
                @test convert(UInt64, RNA_G) == UInt64(2)
                @test convert(UInt64, RNA_U) == UInt64(3)
                @test convert(UInt64, RNA_N) == UInt64(4)
            end
        end

        @testset "Char" begin
            @testset "DNA conversions from Char" begin
                @test convert(DNANucleotide, 'A') == DNA_A
                @test convert(DNANucleotide, 'C') == DNA_C
                @test convert(DNANucleotide, 'G') == DNA_G
                @test convert(DNANucleotide, 'T') == DNA_T
                @test convert(DNANucleotide, 'N') == DNA_N
            end

            @testset "RNA conversions from Char" begin
                @test convert(RNANucleotide, 'A') == RNA_A
                @test convert(RNANucleotide, 'C') == RNA_C
                @test convert(RNANucleotide, 'G') == RNA_G
                @test convert(RNANucleotide, 'U') == RNA_U
                @test convert(RNANucleotide, 'N') == RNA_N
            end

            @testset "DNA conversions to Char" begin
                @test convert(Char, DNA_A) == 'A'
                @test convert(Char, DNA_C) == 'C'
                @test convert(Char, DNA_G) == 'G'
                @test convert(Char, DNA_T) == 'T'
                @test convert(Char, DNA_N) == 'N'
            end

            @testset "RNA conversions to Char" begin
                @test convert(Char, RNA_A) == 'A'
                @test convert(Char, RNA_C) == 'C'
                @test convert(Char, RNA_G) == 'G'
                @test convert(Char, RNA_U) == 'U'
                @test convert(Char, RNA_N) == 'N'
            end
        end

        @testset "Other numeric types" begin
            @test convert(Int, DNA_A) == 0
            @test convert(Int, DNA_C) == 1
            @test convert(Int, DNA_G) == 2
            @test convert(Int, DNA_T) == 3
            @test convert(Int, DNA_N) == 4
            @test convert(DNANucleotide, 0) == DNA_A
            @test convert(DNANucleotide, 1) == DNA_C
            @test convert(DNANucleotide, 2) == DNA_G
            @test convert(DNANucleotide, 3) == DNA_T
            @test convert(DNANucleotide, 4) == DNA_N

            @test convert(Int, RNA_A) == 0
            @test convert(Int, RNA_C) == 1
            @test convert(Int, RNA_G) == 2
            @test convert(Int, RNA_U) == 3
            @test convert(Int, RNA_N) == 4
            @test convert(RNANucleotide, 0) == RNA_A
            @test convert(RNANucleotide, 1) == RNA_C
            @test convert(RNANucleotide, 2) == RNA_G
            @test convert(RNANucleotide, 3) == RNA_U
            @test convert(RNANucleotide, 4) == RNA_N
        end
    end

    @testset "Arithmetic and Order" begin
        @test DNA_A + 1 == DNA_C
        @test DNA_N - 1 == DNA_T
        @test DNA_N - DNA_C == 3
        @test RNA_A + 1 == RNA_C
        @test RNA_N - 1 == RNA_U
        @test RNA_N - RNA_C == 3
        @test DNA_A < DNA_C < DNA_G < DNA_T < DNA_N
        @test RNA_A < RNA_C < RNA_G < RNA_U < RNA_N
        @test !(DNA_A > DNA_G)
        @test !(RNA_A > RNA_G)

        @test collect(alphabet(DNANucleotide)) == [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
        @test collect(alphabet(RNANucleotide)) == [RNA_A, RNA_C, RNA_G, RNA_U, RNA_N]
    end

    @testset "Show DNA" begin
        buf = IOBuffer()
        for nt in [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
            show(buf, nt)
        end
        @test takebuf_string(buf) == "ACGTN"
    end

    @testset "Show RNA" begin
        buf = IOBuffer()
        for nt in [RNA_A, RNA_C, RNA_G, RNA_U, RNA_N]
            show(buf, nt)
        end
        @test takebuf_string(buf) == "ACGUN"
    end

    @testset "Sequences" begin
        function dna_complement(seq::AbstractString)
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

            return convert(AbstractString, seqc)
        end

        function rna_complement(seq::AbstractString)
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

            return convert(AbstractString, seqc)
        end

        function check_reversal(T, seq)
            return reverse(seq) == convert(AbstractString, reverse(convert(T, seq)))
        end

        function check_dna_complement(T, seq)
            return dna_complement(seq) ==
                convert(AbstractString, complement(convert(T, seq)))
        end

        function check_rna_complement(T, seq)
            return rna_complement(seq) ==
                convert(AbstractString, complement(convert(T, seq)))
        end

        function check_dna_revcomp(T, seq)
            return reverse(dna_complement(seq)) ==
                convert(AbstractString, reverse_complement(convert(T, seq)))
        end

        function check_rna_revcomp(T, seq)
            return reverse(rna_complement(seq)) ==
                convert(AbstractString, reverse_complement(convert(T, seq)))
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

        @testset "Nucleotide Sequences" begin
            reps = 10
            @testset "Construction and Conversions" begin
                @testset "Constructing empty sequences" begin
                    # Check construction of empty nucleotide sequences
                    #  using RNASequence and DNASequence functions
                    @test RNASequence() == NucleotideSequence(RNANucleotide)
                    @test DNASequence() == NucleotideSequence(DNANucleotide)
                end

                @testset "Conversion from/to Strings" begin
                    # Check that sequences in strings survive round trip conversion:
                    #   String → NucleotideSequence → String
                    function check_string_construction(T::Type, seq::AbstractString)
                        return convert(AbstractString, NucleotideSequence{T}(seq)) == uppercase(seq)
                    end

                    for len in [0, 1, 10, 32, 1000, 10000, 100000]
                        @test all(Bool[check_string_construction(DNANucleotide, random_dna(len)) for _ in 1:reps])
                        @test all(Bool[check_string_construction(RNANucleotide, random_rna(len)) for _ in 1:reps])
                        @test all(Bool[check_string_construction(DNANucleotide, lowercase(random_dna(len))) for _ in 1:reps])
                        @test all(Bool[check_string_construction(RNANucleotide, lowercase(random_rna(len))) for _ in 1:reps])
                    end

                    # Non-nucleotide characters should throw
                    @test_throws Exception DNASequence("ACCNNCATTTTTTAGATXATAG")
                    @test_throws Exception RNASequence("ACCNNCATTTTTTAGATXATAG")
                end

                @testset "Conversion between RNA and DNA" begin
                    @test convert(RNASequence, DNASequence("ACGTN")) == rna"ACGUN"
                    @test convert(DNASequence, RNASequence("ACGUN")) == dna"ACGTN"
                end

                @testset "Construction from nucleotide vectors" begin
                    function check_vector_construction(T::Type, seq::AbstractString)
                        xs = T[convert(T, c) for c in seq]
                        return NucleotideSequence{T}(xs) == NucleotideSequence{T}(seq)
                    end

                    for len in [0, 1, 10, 32, 1000, 10000, 100000]
                        @test all(Bool[check_vector_construction(DNANucleotide, random_dna(len)) for _ in 1:reps])
                        @test all(Bool[check_vector_construction(RNANucleotide, random_rna(len)) for _ in 1:reps])
                    end
                end

                @testset "Concatenation" begin
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

                        return convert(AbstractString, seq) == uppercase(str)
                    end

                    @test all(Bool[check_concatenation(DNANucleotide, rand(1:10)) for _ in 1:100])
                end

                @testset "Repetition" begin
                    function check_repetition(::Type{DNANucleotide}, n)
                        chunk = random_dna(rand(100:300))
                        start = rand(1:length(chunk))
                        stop = rand(start:length(chunk))

                        str = chunk[start:stop] ^ n
                        seq = DNASequence(chunk)[start:stop] ^ n
                        return convert(AbstractString, seq) == uppercase(str)
                    end

                    @test all(Bool[check_repetition(DNANucleotide, rand(1:10)) for _ in 1:100])
                end
            end

            @testset "Equality" begin
                reps = 10
                function check_seq_equality(len)
                    a = random_dna(len)
                    return a == copy(a)
                end

                for len in [1, 10, 32, 1000]
                    @test all(Bool[check_seq_equality(len) for _ in 1:reps])
                end

                a = b = dna"ACTGN"
                @test a == b
                @test dna"ACTGN" == dna"ACTGN"
                @test dna"ACTGN" != dna"ACTGA"
                @test dna"ACTGN" != dna"ACTG"
                @test dna"ACTG"  != dna"ACTGN"

                c = d = rna"ACUGN"
                @test c == d
                @test rna"ACUGN" == rna"ACUGN"
                @test rna"ACUGN" != rna"ACUGA"
                @test rna"ACUGN" != rna"ACUG"
                @test rna"ACUG"  != rna"ACUGN"

                a = dna"ACGTNACGTN"
                b = dna"""
                ACGTN
                ACGTN
                """
                @test a == b

                c = rna"ACUGNACUGN"
                d = rna"""
                ACUGN
                ACUGN
                """
                @test c == d
            end

            @testset "Length" begin
                for len in [0, 1, 10, 32, 1000, 10000, 100000]
                    @test length(DNASequence(random_dna(len))) == len
                    @test endof(DNASequence(random_dna(len))) == len
                    @test length(DNASequence(random_dna(len))) == len
                    @test endof(DNASequence(random_dna(len))) == len
                end
            end

            @testset "Copy" begin
                function check_copy(T, seq)
                    return convert(AbstractString, copy(NucleotideSequence{T}(seq))) == seq
                end

                for len in [1, 10, 32, 1000, 10000, 100000]
                    @test all(Bool[check_copy(DNANucleotide, random_dna(len)) for _ in 1:reps])
                    @test all(Bool[check_copy(RNANucleotide, random_rna(len)) for _ in 1:reps])
                end
            end

            @testset "Access and Iterations" begin
                dna_seq = dna"ACTG"
                rna_seq = rna"ACUG"

                @testset "Access DNA Sequence" begin
                    # Access indexes out of bounds
                    @test_throws Exception dna_seq[-1]
                    @test_throws Exception dna_seq[0]
                    @test_throws Exception dna_seq[5]
                    @test_throws Exception getindex(dna_seq,-1)
                    @test_throws Exception getindex(dna_seq, 0)
                    @test_throws Exception getindex(dna_seq, 5)
                end

                @testset "Iteration through DNA Sequence" begin
                    @test start(dna"ACNTG") == (1,3)
                    @test start(dna"")      == (1,1)

                    @test next(dna"ACTGN", (2,5)) == (DNA_C, (3,5))
                    @test next(dna"ACTGN", (5,5)) == (DNA_N, (6,6))

                    @test  done(dna"", (1,1))
                    @test !done(dna"ACTGN", (2,5))
                    @test !done(dna"ACTGN", (5,5))
                    @test  done(dna"ACTGN", (6,5))
                    @test !done(dna"ACTGN", (0,5))

                    dna_vector = [DNA_A, DNA_C, DNA_T, DNA_G]
                    @test all(Bool[nucleotide == dna_vector[i] for (i, nucleotide) in enumerate(dna_seq)])
                end

                @testset "Access RNA Sequence" begin
                    # Access indexes out of bounds
                    @test_throws Exception rna_seq[-1]
                    @test_throws Exception rna_seq[0]
                    @test_throws Exception rna_seq[5]
                    @test_throws Exception getindex(rna_seq, -1)
                    @test_throws Exception getindex(rna_seq, 0)
                    @test_throws Exception getindex(rna_seq, 5)
                end

                @testset "Iteration through RNA Sequence" begin
                    @test start(rna"ACNUG") == (1,3)
                    @test start(rna"")      == (1,1)

                    @test next(rna"ACUGN", (2,5)) == (RNA_C, (3,5))
                    @test next(rna"ACUGN", (5,5)) == (RNA_N, (6,6))

                    @test  done(rna"", (1,1))
                    @test !done(rna"ACUGN", (2,5))
                    @test !done(rna"ACUGN", (5,5))
                    @test  done(rna"ACUGN", (6,5))
                    @test !done(rna"ACUGN", (0,5))

                    # Iteration through RNA Sequence

                    rna_vector = [RNA_A, RNA_C, RNA_U, RNA_G]
                    @test all(Bool[nucleotide == rna_vector[i] for (i, nucleotide) in enumerate(rna_seq)])
                end

                @testset "Indexing with Ranges" begin
                    @test getindex(dna"ACTGNACTGN", 1:5) == dna"ACTGN"
                    @test getindex(rna"ACUGNACUGN", 1:5) == rna"ACUGN"
                    @test getindex(dna"ACTGNACTGN", 5:1) == dna""
                    @test getindex(rna"ACUGNACUGN", 5:1) == rna""
                end
            end

            @testset "Subsequence Construction" begin
                for len in [1, 10, 32, 1000, 10000, 100000]
                    seq = random_dna(len)
                    dnaseq = DNASequence(seq)

                    results = Bool[]
                    for _ in 1:reps
                        part = random_interval(1, length(seq))
                        push!(results, seq[part] == convert(AbstractString, dnaseq[part]))
                    end
                    @test all(results)
                end

                for len in [1, 10, 32, 1000, 10000, 100000]
                    seq = random_rna(len)
                    rnaseq = RNASequence(seq)

                    results = Bool[]
                    for _ in 1:reps
                        part = random_interval(1, length(seq))

                        push!(results, seq[part] == convert(AbstractString, rnaseq[part]))
                    end
                    @test all(results)
                end

                @testset "Subsequence Construction from Ranges" begin
                    # Subsequence from range
                    @test RNASequence(rna"AUCGAUCG", 5:8) == RNASequence("AUCG")
                    @test DNASequence(dna"ATCGATCG", 5:8) == DNASequence("ATCG")

                    # Invalid ranges
                    @test_throws Exception RNASequence(rna"AUCGAUCG", 5:10)
                    @test_throws Exception DNASequence(dna"ATCGATCG", 5:10)

                    # Empty ranges
                    @test RNASequence(rna"AUCGAUCG", 5:4) == RNASequence()
                    @test DNASequence(dna"ATCGATCG", 5:4) == DNASequence()

                    # Subsequence of subsequence
                    @test dna"ACGTAG"[4:end][1:2] == dna"TA"
                    @test dna"ACGTAG"[4:end][2:3] == dna"AG"
                    @test_throws Exception dna"ACGTAG"[4:end][0:1]
                    @test_throws Exception dna"ACGTAG"[4:end][3:4]
                end
            end

            @testset "Transformations" begin
                @testset "Reversal" begin
                    for len in [0, 1, 10, 32, 1000, 10000, 100000]
                        @test all(Bool[check_reversal(DNASequence, random_dna(len)) for _ in 1:reps])
                        @test all(Bool[check_reversal(RNASequence, random_rna(len)) for _ in 1:reps])
                    end
                end

                @testset "Complement" begin
                    for len in [1, 10, 32, 1000, 10000, 100000]
                        @test all(Bool[check_dna_complement(DNASequence, random_dna(len)) for _ in 1:reps])
                        @test all(Bool[check_rna_complement(RNASequence, random_rna(len)) for _ in 1:reps])
                    end
                end

                @testset "Reverse Complement" begin
                    for len in [1, 10, 32, 1000, 10000, 100000]
                        @test all(Bool[check_dna_revcomp(DNASequence, random_dna(len)) for _ in 1:reps])
                        @test all(Bool[check_rna_revcomp(RNASequence, random_rna(len)) for _ in 1:reps])
                    end
                end
            end

            @testset "Mismatches" begin
                for len in [1, 10, 32, 1000, 10000, 100000]
                    @test all(Bool[check_mismatches(DNASequence, random_dna(len), random_dna(len)) for _ in 1:reps])
                    @test all(Bool[check_mismatches(RNASequence, random_rna(len), random_rna(len)) for _ in 1:reps])
                end
            end

            @testset "Mutability" begin
                seq = dna"ACGTACGT"
                @test !ismutable(seq)
                @test_throws Exception seq[1] = DNA_C

                seq2 = seq[1:4]
                @test !ismutable(seq2)

                mutable!(seq)
                @test ismutable(seq)

                seq[1] = DNA_C
                @test seq[1] == DNA_C
                @test seq2[1] == DNA_A

                seq[2] = DNA_N
                @test seq[2] == DNA_N
                @test seq2[2] == DNA_C

                seq[2] = DNA_G
                @test seq[2] == DNA_G

                immutable!(seq)
                @test_throws Exception seq[1] = DNA_A

                seq = dna"ACGTACGT"
                mutable!(seq)
                rnaseq = convert(RNASequence, seq)
                @test ismutable(rnaseq)
                @test_throws Exception rnaseq[1] = DNA_C
                rnaseq[1] = RNA_C
                @test seq == dna"ACGTACGT"
                @test rnaseq == rna"CCGUACGU"
            end
        end

        @testset "SequenceNIterator" begin
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
                @test all(Bool[check_ns(DNASequence, random_dna(len)) for _ in 1:reps])
                @test all(Bool[check_ns(RNASequence, random_rna(len)) for _ in 1:reps])
            end

            dna_seq   = dna"ANANANA"
            dna_niter = npositions(dna_seq)
            @test start(dna_niter) == 2
            @test next(dna_niter, 2) == (2,4)
            @test !done(dna_niter, 6)
            @test  done(dna_niter, 8)


            ns = [2,4,6]
            for (i, n) in enumerate(dna_niter)
                @test n == ns[i]
            end

            rna_seq   = rna"ANANANA"
            rna_niter = npositions(rna_seq)
            @test start(rna_niter) == 2
            @test next(rna_niter, 2) == (2,4)
            @test !done(rna_niter, 6)
            @test  done(rna_niter, 8)

            ns = [2,4,6]
            for (i, n) in enumerate(rna_niter)
                @test n == ns[i]
            end

        end

        @testset "Kmer" begin
            reps = 100
            @testset "Construction and Conversions" begin
                # Check that kmers in strings survive round trip conversion:
                #   UInt64 → Kmer → UInt64
                function check_uint64_convertion(T::Type, n::UInt64, len::Int)
                    return convert(UInt64, convert(Kmer{T, len}, n)) === n
                end

                # Check that kmers in strings survive round trip conversion:
                #   String → Kmer → String
                function check_string_construction(T::Type, seq::AbstractString)
                    return convert(AbstractString, convert(Kmer{T}, seq)) == uppercase(seq)
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
                    return convert(AbstractString, [convert(Char, c) for c in seq]) ==
                           convert(AbstractString, kmer(seq...))
                end

                # Check that kmers in strings survive round trip conversion:
                #   String → NucleotideSequence → Kmer → NucleotideSequence → String
                function check_roundabout_construction(T::Type, seq::AbstractString)
                    return convert(AbstractString,
                               convert(NucleotideSequence{T},
                                   convert(Kmer,
                                       convert(NucleotideSequence{T}, seq)))) == uppercase(seq)
                end



                for len in [0, 1, 16, 32]
                    # UInt64 conversions
                    @test all(Bool[check_uint64_convertion(DNANucleotide, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])
                    @test all(Bool[check_uint64_convertion(RNANucleotide, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])

                    # String construction
                    @test all(Bool[check_string_construction(DNANucleotide, random_dna_kmer(len)) for _ in 1:reps])
                    @test all(Bool[check_string_construction(RNANucleotide, random_rna_kmer(len)) for _ in 1:reps])

                    # DNA/RNASequence Constructions
                    @test all(Bool[check_dnasequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps])
                    @test all(Bool[check_rnasequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps])

                    # NucleotideSequence Construction
                    @test all(Bool[check_nucsequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps])
                    @test all(Bool[check_nucsequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps])

                    # Construction from nucleotide arrays
                    if len > 0
                        @test all(Bool[check_nucarray_kmer(random_dna_kmer_nucleotides(len)) for _ in 1:reps])
                        @test all(Bool[check_nucarray_kmer(random_rna_kmer_nucleotides(len)) for _ in 1:reps])
                    end

                    # Roundabout conversions
                    @test all(Bool[check_roundabout_construction(DNANucleotide, random_dna_kmer(len)) for _ in 1:reps])
                    @test all(Bool[check_roundabout_construction(RNANucleotide, random_rna_kmer(len)) for _ in 1:reps])
                end

                @test_throws Exception kmer() # can't construct 0-mer using `kmer()`
                @test_throws Exception kmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
                @test_throws Exception kmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
                @test_throws Exception kmer(rna"ACGNU")# no Ns in kmers
                @test_throws Exception rnmakmer(rna"ACGNU")# no Ns in kmers
                @test_throws Exception kmer(dna"ACGNT") # no Ns in kmers
                @test_throws Exception dnmakmer(dna"ACGNT") # no Ns in kmers
                @test_throws Exception kmer(RNA_A, DNA_A) # no mixing of RNA and DNA
                @test_throws Exception kmer(random_rna(33)) # no kmer larger than 32nt
                @test_throws Exception kmer(random_dna(33)) # no kmer larger than 32nt
                @test_throws Exception kmer(RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 32nt
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U,
                                  RNA_A, RNA_C, RNA_G, RNA_U)
                @test_throws Exception kmer(DNA_A, DNA_C, DNA_G, DNA_T, # no kmer larger than 32nt
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T,
                                  DNA_A, DNA_C, DNA_G, DNA_T)

                @testset "From strings" begin
                    @test dnakmer("ACTG") == convert(Kmer, DNASequence("ACTG"))
                    @test rnakmer("ACUG") == convert(Kmer, RNASequence("ACUG"))

                    # N is not allowed in Kmers
                    @test_throws Exception dnakmer("ACGTNACGT")
                    @test_throws Exception rnakmer("ACGUNACGU")
                end
            end

            @testset "Comparisons" begin
                @testset "Equality" begin
                    function check_seq_kmer_equality(len)
                        a = dnakmer(random_dna_kmer(len))
                        b = convert(DNASequence, a)
                        return a == b && b == a
                    end

                    for len in [1, 10, 32]
                        @test all(Bool[check_seq_kmer_equality(len) for _ in 1:reps])
                    end

                    # True negatives
                    @test (dnakmer("ACG") == rnakmer("ACG")) == false
                    @test (dnakmer("T")   == rnakmer("U"))   == false
                    @test (dnakmer("AC")  == dnakmer("AG"))  == false
                    @test (rnakmer("AC")  == rnakmer("AG"))  == false

                    @test (dnakmer("ACG") == rna"ACG") == false
                    @test (dnakmer("T")   == rna"U")   == false
                    @test (dnakmer("AC")  == dna"AG")  == false
                    @test (rnakmer("AC")  == rna"AG")  == false

                    @test (rna"ACG" == dnakmer("ACG")) == false
                    @test (rna"U"   == dnakmer("T"))   == false
                    @test (dna"AG"  == dnakmer("AC"))  == false
                    @test (rna"AG"  == rnakmer("AC"))  == false
                end

                @testset "Inequality" begin
                    for len in [1, 10, 32]
                        @test  isless(convert(DNAKmer{1}, UInt64(0)), convert(DNAKmer{1}, UInt64(1)))
                        @test !isless(convert(DNAKmer{1}, UInt64(0)), convert(DNAKmer{1}, UInt64(0)))
                        @test !isless(convert(DNAKmer{1}, UInt64(1)), convert(DNAKmer{1}, UInt64(0)))

                        @test  isless(convert(RNAKmer{1}, UInt64(0)), convert(RNAKmer{1}, UInt64(1)))
                        @test !isless(convert(RNAKmer{1}, UInt64(0)), convert(RNAKmer{1}, UInt64(0)))
                        @test !isless(convert(RNAKmer{1}, UInt64(1)), convert(RNAKmer{1}, UInt64(0)))
                    end
                end

                @testset "Hash" begin
                    kmers = map(dnakmer, ["AAAA", "AACT", "ACGT", "TGCA"])
                    for x in kmers, y in kmers
                        @test (x == y) == (hash(x) == hash(y))
                    end
                    kmers = map(rnakmer, ["AAAA", "AACU", "ACGU", "UGCA"])
                    for x in kmers, y in kmers
                        @test (x == y) == (hash(x) == hash(y))
                    end
                end
            end

            @testset "Length" begin
                for len in [0, 1, 16, 32]
                    @test length(dnakmer(random_dna_kmer(len))) == len
                    @test length(rnakmer(random_rna_kmer(len))) == len
                end
            end

            @testset "Arithmetic" begin
                x = dnakmer("AA")
                @test x - 1 == x + (-1) == dnakmer("TT")
                @test x + 1 == x - (-1) == dnakmer("AC")

                x = dnakmer("TT")
                @test x - 1 == x + (-1) == dnakmer("TG")
                @test x + 1 == x - (-1) == dnakmer("AA")

                base = dnakmer("AAA")
                offset = 0
                nucs = "ACGT"
                for a in nucs, b in nucs, c in nucs
                    @test base + offset == dnakmer(string(a, b, c))
                    offset += 1
                end
            end

            @testset "Order" begin
                @test dnakmer("AA") < dnakmer("AC") < dnakmer("AG") < dnakmer("AT") < dnakmer("CA")
                @test rnakmer("AA") < rnakmer("AC") < rnakmer("AG") < rnakmer("AU") < rnakmer("CA")
            end

            @testset "Access and Iterations" begin
                dna_kmer = dnakmer("ACTG")
                rna_kmer = rnakmer("ACUG")

                @testset "Access DNA Kmer" begin
                    @test dna_kmer[1] == DNA_A
                    @test dna_kmer[2] == DNA_C
                    @test dna_kmer[3] == DNA_T
                    @test dna_kmer[4] == DNA_G

                    # Access indexes out of bounds
                    @test_throws Exception dna_kmer[-1]
                    @test_throws Exception dna_kmer[0]
                    @test_throws Exception dna_kmer[5]
                    @test_throws Exception getindex(dna_kmer,-1)
                    @test_throws Exception getindex(dna_kmer, 0)
                    @test_throws Exception getindex(dna_kmer, 5)
                end

                @testset "Iteration through DNA Kmer" begin
                    @test start(dnakmer("ACTG"))  == 1
                    @test start(dnakmer(""))      == 1

                    @test next(dnakmer("ACTG"), 1) == (DNA_A, 2)
                    @test next(dnakmer("ACTG"), 4) == (DNA_G, 5)

                    @test  done(dnakmer(""), 1)
                    @test !done(dnakmer("ACTG"), 1)
                    @test !done(dnakmer("ACTG"), 4)
                    @test  done(dnakmer("ACTG"), 5)
                    @test !done(dnakmer("ACTG"), -1)


                    dna_kmer_vector = [DNA_A, DNA_C, DNA_T, DNA_G]
                    @test all(Bool[nucleotide == dna_kmer_vector[i] for (i, nucleotide) in enumerate(dna_kmer)])
                end

                @testset "Access RNA Kmer" begin
                    @test rna_kmer[1] == RNA_A
                    @test rna_kmer[2] == RNA_C
                    @test rna_kmer[3] == RNA_U
                    @test rna_kmer[4] == RNA_G

                    # Access indexes out of bounds
                    @test_throws Exception rna_kmer[-1]
                    @test_throws Exception rna_kmer[0]
                    @test_throws Exception rna_kmer[5]
                    @test_throws Exception getindex(rna_kmer, -1)
                    @test_throws Exception getindex(rna_kmer, 0)
                    @test_throws Exception getindex(rna_kmer, 5)
                end

                @testset "Iteration through RNA Kmer" begin
                    @test start(rnakmer("ACUG"))  == 1
                    @test start(rnakmer(""))      == 1

                    @test next(rnakmer("ACUG"), 1) == (RNA_A, 2)
                    @test next(rnakmer("ACUG"), 4) == (RNA_G, 5)

                    @test  done(rnakmer(""), 1)
                    @test !done(rnakmer("ACUG"), 1)
                    @test !done(rnakmer("ACUG"), 4)
                    @test  done(rnakmer("ACUG"), 5)
                    @test !done(rnakmer("ACUG"), -1)

                    rna_kmer_vector = [RNA_A, RNA_C, RNA_U, RNA_G]
                    @test all(Bool[nucleotide == rna_kmer_vector[i] for (i, nucleotide) in enumerate(rna_kmer)])
                end
            end

            reps = 10
            @testset "Transformations" begin
                @testset "Reversal" begin
                    for len in [0, 1, 16, 32]
                        @test all(Bool[check_reversal(DNAKmer, random_dna_kmer(len)) for _ in 1:reps])
                        @test all(Bool[check_reversal(RNAKmer, random_rna_kmer(len)) for _ in 1:reps])
                    end
                end

                @testset "Complement" begin
                    for len in [0, 1, 16, 32]
                        @test all(Bool[check_dna_complement(DNAKmer, random_dna_kmer(len)) for _ in 1:reps])
                        @test all(Bool[check_rna_complement(RNAKmer, random_rna_kmer(len)) for _ in 1:reps])
                    end
                end

                @testset "Reverse Complement" begin
                    for len in [0, 1, 16, 32]
                        @test all(Bool[check_dna_revcomp(DNAKmer, random_dna_kmer(len)) for _ in 1:reps])
                        @test all(Bool[check_rna_revcomp(RNAKmer, random_rna_kmer(len)) for _ in 1:reps])
                    end
                end
            end

            @testset "Mismatches" begin
                for len in [0, 1, 16, 32]
                    @test all(Bool[check_mismatches(DNAKmer, random_dna_kmer(len), random_dna_kmer(len)) for _ in 1:reps])
                    @test all(Bool[check_mismatches(RNAKmer, random_rna_kmer(len), random_rna_kmer(len)) for _ in 1:reps])
                end
            end
        end

        @testset "EachKmer" begin
            function string_eachkmer(seq::AbstractString, k, step)
                kmers = AbstractString[]
                i = 1
                for i in 1:step:length(seq) - k + 1
                    subseq = seq[i:i + k - 1]
                    if !in('N', subseq)
                        push!(kmers, subseq)
                    end
                end
                return kmers
            end

            function check_eachkmer(T, seq::AbstractString, k, step)
                xs = [convert(AbstractString, x) for (i, x) in collect(each(Kmer{T, k}, NucleotideSequence{T}(seq), step))]
                ys = [convert(AbstractString, x) for (i, x) in collect(eachkmer(NucleotideSequence{T}(seq), k, step))]
                zs = string_eachkmer(seq, k, step)
                return xs == ys == zs
            end

            reps = 10
            len = 10000

            for k in [0, 1, 3, 16, 32], step in 1:3
                @test all(Bool[check_eachkmer(DNANucleotide, random_dna(len), k, step) for _ in 1:reps])
                @test all(Bool[check_eachkmer(RNANucleotide, random_rna(len), k, step) for _ in 1:reps])
            end

            @test isempty(collect(each(DNAKmer{1}, dna"")))
            @test isempty(collect(each(DNAKmer{1}, dna"NNNNNNNNNN")))
            @test_throws Exception each(DNAKmer{-1}, dna"ACGT")
            @test_throws Exception each(DNAKmer{33}, dna"ACGT")
        end

        @testset "Nucleotide Counting" begin
            function string_nucleotide_count(::Type{DNANucleotide}, seq::AbstractString)
                counts = Dict{DNANucleotide, Int}(
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

            function string_nucleotide_count(::Type{RNANucleotide}, seq::AbstractString)
                counts = Dict{RNANucleotide, Int}(
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

            function check_nucleotide_count(::Type{DNANucleotide}, seq::AbstractString)
                string_counts = string_nucleotide_count(DNANucleotide, seq)
                seq_counts = NucleotideCounts(DNASequence(seq))
                return string_counts[DNA_A] == seq_counts[DNA_A] &&
                       string_counts[DNA_C] == seq_counts[DNA_C] &&
                       string_counts[DNA_G] == seq_counts[DNA_G] &&
                       string_counts[DNA_T] == seq_counts[DNA_T] &&
                       string_counts[DNA_N] == seq_counts[DNA_N]
            end

            function check_nucleotide_count(::Type{RNANucleotide}, seq::AbstractString)
                string_counts = string_nucleotide_count(RNANucleotide, seq)
                seq_counts = NucleotideCounts(RNASequence(seq))
                return string_counts[RNA_A] == seq_counts[RNA_A] &&
                       string_counts[RNA_C] == seq_counts[RNA_C] &&
                       string_counts[RNA_G] == seq_counts[RNA_G] &&
                       string_counts[RNA_U] == seq_counts[RNA_U] &&
                       string_counts[RNA_N] == seq_counts[RNA_N]
            end

            function check_kmer_nucleotide_count(::Type{DNANucleotide}, seq::AbstractString)
                string_counts = string_nucleotide_count(DNANucleotide, seq)
                kmer_counts = NucleotideCounts(dnakmer(seq))
                return string_counts[DNA_A] == kmer_counts[DNA_A] &&
                       string_counts[DNA_C] == kmer_counts[DNA_C] &&
                       string_counts[DNA_G] == kmer_counts[DNA_G] &&
                       string_counts[DNA_T] == kmer_counts[DNA_T] &&
                       string_counts[DNA_N] == kmer_counts[DNA_N]
            end

            function check_kmer_nucleotide_count(::Type{RNANucleotide}, seq::AbstractString)
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
                @test all(Bool[check_nucleotide_count(DNANucleotide, random_dna(len)) for _ in 1:reps])
                @test all(Bool[check_nucleotide_count(RNANucleotide, random_rna(len)) for _ in 1:reps])
            end

            for len in [1, 10, 32]
                @test all(Bool[check_kmer_nucleotide_count(DNANucleotide, random_dna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_kmer_nucleotide_count(RNANucleotide, random_rna_kmer(len)) for _ in 1:reps])
            end
        end

        @testset "Kmer Counting" begin
            function string_kmer_count{T <: Nucleotide}(::Type{T}, seq::AbstractString, k, step)
                counts = Dict{Kmer{T, k}, Int}()
                for x in UInt64(0):UInt64(4^k-1)
                    counts[convert(Kmer{T, k}, x)] = 0
                end

                for i in 1:step:(length(seq)-k+1)
                    s = seq[i:i+k-1]
                    if 'N' in s
                        continue
                    end
                    counts[convert(Kmer{T, k}, s)] += 1
                end

                return counts
            end

            function check_kmer_count{T <: Nucleotide}(::Type{T}, seq::AbstractString, k, step)
                string_counts = string_kmer_count(T, seq, k, step)
                kmer_counts = KmerCounts{T, k}(convert(NucleotideSequence{T}, seq), step)
                for y in UInt64(0):UInt64(4^k-1)
                    x = convert(Kmer{T, k}, y)
                    if string_counts[x] != kmer_counts[x]
                        return false
                    end
                end
                return true
            end

            reps = 10
            for len in [1, 10, 32, 1000, 10000]
                for k in [1, 2, 5]
                    for step in [1, 3]
                        @test all(Bool[check_kmer_count(DNANucleotide, random_dna(len), k, step) for _ in 1:reps])
                        @test all(Bool[check_kmer_count(RNANucleotide, random_rna(len), k, step) for _ in 1:reps])
                    end
                end
            end
        end
    end
end

@testset "Aminoacids" begin
    @testset "Arithmetic and Order" begin
        @test AA_A + 1 == AA_R
        @test AA_R + 1 == AA_N
        @test AA_R - 1 == AA_A
        @test AA_D - AA_A ==  3
        @test AA_A - AA_D == -3
        @test AA_A < AA_R < AA_N < AA_Z < AA_X < AA_O < AA_U
        @test !(AA_J < AA_B)

        @test length(alphabet(AminoAcid)) == 26
        @test AA_A in alphabet(AminoAcid)
        @test AA_I in alphabet(AminoAcid)
        @test AA_U in alphabet(AminoAcid)
    end

    @testset "AminoAcid Sequences" begin
        reps = 10
        @testset "Construction" begin
            # Non-aa characters should throw
            @test_throws Exception AminoAcidSequence("ATGHLMY@ZACAGNM")

            # Check that sequences in strings survive round trip conversion:
            #   String → AminoAcidSequence → String
            function check_string_construction(seq::AbstractString)
                return convert(AbstractString, AminoAcidSequence(seq)) == uppercase(seq)
            end

            for len in [0, 1, 10, 32, 1000, 10000, 100000]
                @test all(Bool[check_string_construction(random_aa(len)) for _ in 1:reps])
                @test all(Bool[check_string_construction(lowercase(random_aa(len))) for _ in 1:reps])
            end

            # Check creation of empty
        end

        @testset "Conversion" begin
            seq = aa"ARNDCQEGHILKMFPSTWYVOUBZJX"
            @test convert(AminoAcidSequence, [aa for aa in seq]) == seq
            @test convert(Vector{AminoAcid}, seq) == [aa for aa in seq]

            @test_throws Exception convert(AminoAcidSequence, [convert(AminoAcid, UInt8(30))])
        end

        @testset "Concatenation" begin
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

                return convert(AbstractString, seq) == uppercase(str)
            end

            @test all(Bool[check_concatenation(rand(1:10)) for _ in 1:100])
        end

        @testset "Equality" begin
            seq = aa"ARNDCQEGHILKMFPSTWYVX"
            @test seq == aa"ARNDCQEGHILKMFPSTWYVX"
            @test seq != aa"ARNDCQEGHILKMFPSTWYXV"
            @test seq != aa"ARNDCQEGHLKMFPSTWYVX"

            seq′ = aa"""
            ARNDCQEGHI
            LKMFPSTWYV
            X
            """
            @test seq == seq′
        end

        @testset "Repetition" begin
            function check_repetition(n)
                chunk = random_aa(rand(100:300))
                start = rand(1:length(chunk))
                stop = rand(start:length(chunk))

                str = chunk[start:stop] ^ n
                seq = AminoAcidSequence(chunk)[start:stop] ^ n

                return convert(AbstractString, seq) == uppercase(str)
            end

            @test all(Bool[check_repetition(rand(1:10)) for _ in 1:100])
        end

        @testset "Copy" begin
            function check_copy(seq)
                return convert(AbstractString, copy(AminoAcidSequence(seq))) == uppercase(seq)
            end

            for len in [1, 10, 32, 1000, 10000, 100000]
                @test all(Bool[check_copy(random_aa(len)) for _ in 1:reps])
            end
        end

        @testset "Subsequence Construction" begin
            for len in [1, 10, 32, 1000, 10000, 100000]
                seq = random_aa(len)
                aaseq = AminoAcidSequence(seq)

                results = Bool[]
                for _ in 1:reps
                    part = random_interval(1, length(seq))
                    push!(results, seq[part] == convert(AbstractString, aaseq[part]))
                end
                @test all(results)
            end
        end

        @testset "Mutability" begin
            seq = aa"ARNDCQEGHILKMFPSTWYVX"
            @test ismutable(seq) == false
            @test_throws Exception seq[1] = AA_C

            seq2 = seq[1:4]
            @test ismutable(seq2) == false

            mutable!(seq)
            @test ismutable(seq)
            seq[1] = AA_C
            @test seq[1] == AA_C
            @test seq2[1] == AA_A

            immutable!(seq)
            @test_throws Exception seq[1] = AA_A
        end
    end

    @testset "Parsers" begin
        @testset "Valid Cases" begin
            # case-insensitive and ignores spaces
            @test parse(AminoAcid, "a") == AA_A
            @test parse(AminoAcid, "Ala") == AA_A
            @test parse(AminoAcid, "aLa ") == AA_A
            @test parse(AminoAcid, " alA ") == AA_A
            @test parse(AminoAcid, "\tAlA\n") == AA_A
            @test parse(AminoAcid, "x") == AA_X
            @test parse(AminoAcid, "X") == AA_X
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
                ("B", "ASX", AA_B),
                ("J", "XLE", AA_J),
                ("Z", "GLX", AA_Z),
                ("X", "XAA", AA_X),
                ("O", "PYL", AA_O),
                ("U", "SEC", AA_U),
            ]
            @test length(aas) == 26
            for (one, three, aa) in aas
                @test parse(AminoAcid, one) == aa
                @test parse(AminoAcid, three) == aa
            end
        end

        @testset "Invalid Cases" begin
            @test_throws Exception parse(AminoAcid, "")
            @test_throws Exception parse(AminoAcid, "AL")
            @test_throws Exception parse(AminoAcid, "LA")
            @test_throws Exception parse(AminoAcid, "ALAA")
        end
    end
end

@testset "Translation" begin
    # crummy string translation to test against
    standard_genetic_code_dict = Dict{AbstractString, Char}(
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

    function string_translate(seq::AbstractString)
        @assert length(seq) % 3 == 0
        aaseq = Array(Char, div(length(seq), 3))
        for i in 1:3:length(seq) - 3 + 1
            aaseq[div(i, 3) + 1] = standard_genetic_code_dict[seq[i:i+2]]
        end
        return convert(AbstractString, aaseq)
    end

    function check_translate(seq::AbstractString)
        return string_translate(seq) == convert(AbstractString, translate(RNASequence(seq)))
    end

    reps = 10
    #for len in [1, 10, 32, 1000, 10000, 100000]
    for len in [1, 10, 32]
        @test all(Bool[check_translate(random_translatable_rna(len)) for _ in 1:reps])
    end

    @test_throws Exception translate(dna"ACGTACGTA") # can't translate DNA
    @test_throws Exception translate(rna"ACGUACGU")  # can't translate non-multiples of three
    @test_throws Exception translate(rna"ACGUACGNU") # can't translate N
end


@testset "Sequence Parsing" begin
    @testset "FASTA Parsing" begin
        get_bio_fmt_specimens()

        function check_fasta_parse(filename)
            # Reading from a stream
            for seqrec in open(filename, FASTA)
            end

            # Reading from a memory mapped file
            for seqrec in open(filename, FASTA, memory_map=true)
            end

            # Check round trip
            output = IOBuffer()
            expected_entries = Any[]
            for seqrec in open(filename, FASTA)
                write(output, seqrec)
                push!(expected_entries, seqrec)
            end

            read_entries = Any[]
            for seqrec in open(takebuf_array(output), FASTA)
                push!(read_entries, seqrec)
            end

            return expected_entries == read_entries
        end

        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "FASTA")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            tags = specimen["tags"]
            valid = get(specimen, "valid", true)
            # currently unsupported features
            if contains(tags, "gaps") || contains(tags, "comments") || (contains(tags, "ambiguity") && !contains(tags, "protein"))
                continue
            end
            if valid
                @test check_fasta_parse(joinpath(path, specimen["filename"]))
            else
                @test_throws Exception check_fasta_parse(joinpath(path, specimen["filename"]))
            end
        end
    end

    @testset "FASTQ Parsing" begin
        get_bio_fmt_specimens()

        function check_fastq_parse(filename)
            # Reading from a stream
            for seqrec in open(filename, FASTQ)
            end

            # Reading from a memory mapped file
            for seqrec in open(filename, FASTQ, memory_map=true)
            end

            # Check round trip
            output = IOBuffer()
            expected_entries = Any[]
            for seqrec in open(filename, FASTQ)
                write(output, seqrec)
                push!(expected_entries, seqrec)
            end

            read_entries = Any[]
            for seqrec in open(takebuf_array(output), FASTQ)
                push!(read_entries, seqrec)
            end

            return expected_entries == read_entries
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
                @test check_fastq_parse(joinpath(path, specimen["filename"]))
            else
                @test_throws Exception check_fastq_parse(joinpath(path, specimen["filename"]))
            end
        end
    end
end

@testset "Quality scores" begin
    @testset "Decoding PHRED scores" begin

        function test_decode(encoding, values, expected)
            result = Array(Int8, length(expected))
            Seq.decode_quality_string!(encoding, values, result, 1, length(result))
            @test result == expected

            # Without start and end does the entire string
            Seq.decode_quality_string!(encoding, values, result)
            @test result == expected

            # Non in-place version of the above
            result = Seq.decode_quality_string(encoding, values, 1, length(result))
            @test result == expected
            resutlt = Seq.decode_quality_string(encoding, values)
            @test result == expected
        end

        test_decode(Seq.SANGER_QUAL_ENCODING,
                    UInt8['!', '#', '$', '%', '&', 'I', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 93])

        test_decode(Seq.SOLEXA_QUAL_ENCODING,
                    UInt8[';', 'B', 'C', 'D', 'E', 'h', '~'],
                    Int8[-5, 2, 3, 4, 5, 40, 62])

        test_decode(Seq.ILLUMINA13_QUAL_ENCODING,
                    UInt8['@', 'B', 'C', 'D', 'E', 'h', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 62])

        test_decode(Seq.ILLUMINA15_QUAL_ENCODING,
                    UInt8['C', 'D', 'E', 'h', '~'],
                    Int8[3, 4, 5, 40, 62])

        test_decode(Seq.ILLUMINA18_QUAL_ENCODING,
                    UInt8['!', '#', '$', '%', '&', 'I', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 93])
    end

    @testset "Encoding PHRED scores" begin

        function test_encode(encoding, values, expected)
            # With start & end
            result = Array(UInt8, length(expected))
            Seq.encode_quality_string!(encoding, values, result, 1, length(result))
            @test result == expected
            # Without start & end means the entire length
            Seq.encode_quality_string!(encoding, values, result)
            @test result == expected

            result = Seq.encode_quality_string(encoding, values, 1, length(result))
            @test result == expected
            result = Seq.encode_quality_string(encoding, values)
            @test result == expected
        end

        test_encode(Seq.SANGER_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 93],
                    UInt8['!', '#', '$', '%', '&', 'I', '~'])

        test_encode(Seq.SOLEXA_QUAL_ENCODING,
                    Int8[-5, 2, 3, 4, 5, 40, 62],
                    UInt8[';', 'B', 'C', 'D', 'E', 'h', '~'])

        test_encode(Seq.ILLUMINA13_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 62],
                    UInt8['@', 'B', 'C', 'D', 'E', 'h', '~'])

        test_encode(Seq.ILLUMINA15_QUAL_ENCODING,
                    Int8[3, 4, 5, 40, 62],
                    UInt8['C', 'D', 'E', 'h', '~'])

        test_encode(Seq.ILLUMINA18_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 93],
                    UInt8['!', '#', '$', '%', '&', 'I', '~'])
    end

    @testset "Sequence Writing" begin
        @testset "FASTA writing" begin
            dna_seq1 = random_dna(79 * 3) # full lines
            dna_seq2 = random_dna(50)     # short line
            seq_name1 = "Sequence 1"
            seq_name2 = "Sequence 2"
            seq_description1 = "Description 1"
            seq_description2 = "Description 2"

            wrapped_seq1 = join([dna_seq1[1:79], dna_seq1[80:158], dna_seq1[159:end]], "\n")

            expected_seq1 = string(">", seq_name1, " ", seq_description1, "\n", wrapped_seq1, "\n")
            expected_seq2 = string(">", seq_name2, " ", seq_description2, "\n", dna_seq2, "\n")
            expected = string(expected_seq1, expected_seq2)

            fasta_seq1 = Seq.FASTADNASeqRecord(
                seq_name1, DNASequence(dna_seq1), Seq.FASTAMetadata(seq_description1))
            fasta_seq2 = Seq.FASTADNASeqRecord(
                seq_name2, DNASequence(dna_seq2), Seq.FASTAMetadata(seq_description2))
            sequences = [fasta_seq1, fasta_seq2]

            output = IOBuffer()
            for seq in sequences
                write(output, seq)
            end
            @test takebuf_string(output) == expected
        end
    end
end

end # TestSeq
