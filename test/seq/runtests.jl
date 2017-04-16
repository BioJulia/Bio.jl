module TestSeq

using Base.Test

using Bio: Seq, Var
import BioSymbols
using IntervalTrees.IntervalValue
using BufferedStreams
using StatsBase
using YAML
using TestFunctions
using PairwiseListMatrices

typealias PWM PairwiseListMatrix

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
        "UAA", "UAC", "UAG", "UAU",
        "UCA", "UCC", "UCG", "UCU",
        "UGA", "UGC", "UGG", "UGU",
        "UUA", "UUC", "UUG", "UUU",
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
            seqc[i] = seq[i]
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
            seqc[i] = seq[i]
        end
    end
    return convert(AbstractString, seqc)
end

function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
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

# NOTE: Most tests related to biological symbols are located in BioSymbols.jl.
@testset "Symbols" begin
    @testset "DNA" begin
        @test DNA_A === BioSymbols.DNA_A
        @test ACGT  === BioSymbols.ACGT
        @test ACGTN === BioSymbols.ACGTN
        @test typeof(DNA_A) === BioSymbols.DNA
    end

    @testset "RNA" begin
        @test RNA_A === BioSymbols.RNA_A
        @test ACGU  === BioSymbols.ACGU
        @test ACGUN === BioSymbols.ACGUN
        @test typeof(RNA_A) === BioSymbols.RNA
    end

    @testset "AminoAcid" begin
        @test AA_A === BioSymbols.AA_A
        @test typeof(AA_A) === BioSymbols.AminoAcid
    end

    @testset "Predicate functions" begin
        @test iscompatible(DNA_A, DNA_N)
        @test isambiguous(DNA_N)
        @test iscertain(DNA_A)
        @test isgap(DNA_Gap)
        @test ispurine(DNA_A)
        @test ispyrimidine(DNA_C)
        @test isGC(DNA_G)
    end

    @testset "Misc. functions" begin
        @test length(alphabet(DNA)) == 16
        @test gap(DNA) === DNA_Gap
        @test complement(DNA_A) === DNA_T
    end

    @testset "Encoder" begin
        encode = Seq.encode
        EncodeError = Seq.EncodeError

        @testset "DNA" begin
            # 2 bits
            @test encode(DNAAlphabet{2}, DNA_A) === 0x00
            @test encode(DNAAlphabet{2}, DNA_C) === 0x01
            @test encode(DNAAlphabet{2}, DNA_G) === 0x02
            @test encode(DNAAlphabet{2}, DNA_T) === 0x03
            @test_throws EncodeError encode(DNAAlphabet{2}, DNA_M)
            @test_throws EncodeError encode(DNAAlphabet{2}, DNA_N)
            @test_throws EncodeError encode(DNAAlphabet{2}, DNA_Gap)

            # 4 bits
            for nt in alphabet(DNA)
                @test encode(DNAAlphabet{4}, nt) === reinterpret(UInt8, nt)
            end
            @test_throws EncodeError encode(DNAAlphabet{4}, reinterpret(DNA, 0b10000))
        end

        @testset "RNA" begin
            # 2 bits
            @test encode(RNAAlphabet{2}, RNA_A) === 0x00
            @test encode(RNAAlphabet{2}, RNA_C) === 0x01
            @test encode(RNAAlphabet{2}, RNA_G) === 0x02
            @test encode(RNAAlphabet{2}, RNA_U) === 0x03
            @test_throws EncodeError encode(RNAAlphabet{2}, RNA_M)
            @test_throws EncodeError encode(RNAAlphabet{2}, RNA_N)
            @test_throws EncodeError encode(RNAAlphabet{2}, RNA_Gap)

            # 4 bits
            for nt in alphabet(RNA)
                @test encode(RNAAlphabet{4}, nt) === reinterpret(UInt8, nt)
            end
            @test_throws EncodeError encode(RNAAlphabet{4}, reinterpret(RNA, 0b10000))
        end

        @testset "AminoAcid" begin
            @test encode(AminoAcidAlphabet, AA_A) === 0x00
            for aa in alphabet(AminoAcid)
                @test encode(AminoAcidAlphabet, aa) === convert(UInt8, aa)
            end
            @test_throws Seq.EncodeError encode(AminoAcidAlphabet, BioSymbols.AA_INVALID)
        end
    end

    @testset "Decoder" begin
        decode = Seq.decode
        DecodeError = Seq.DecodeError

        @testset "DNA" begin
            # 2 bits
            @test decode(DNAAlphabet{2}, 0x00) === DNA_A
            @test decode(DNAAlphabet{2}, 0x01) === DNA_C
            @test decode(DNAAlphabet{2}, 0x02) === DNA_G
            @test decode(DNAAlphabet{2}, 0x03) === DNA_T
            @test_throws DecodeError decode(DNAAlphabet{2}, 0x04)
            @test_throws DecodeError decode(DNAAlphabet{2}, 0x0e)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(DNAAlphabet{4}, x) === reinterpret(DNA, x)
            end
            @test_throws DecodeError decode(DNAAlphabet{4}, 0b10000)
        end

        @testset "RNA" begin
            # 2 bits
            @test decode(RNAAlphabet{2}, 0x00) === RNA_A
            @test decode(RNAAlphabet{2}, 0x01) === RNA_C
            @test decode(RNAAlphabet{2}, 0x02) === RNA_G
            @test decode(RNAAlphabet{2}, 0x03) === RNA_U
            @test_throws DecodeError decode(RNAAlphabet{2}, 0x04)
            @test_throws DecodeError decode(RNAAlphabet{2}, 0x0e)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(RNAAlphabet{4}, x) === reinterpret(RNA, x)
            end
            @test_throws DecodeError decode(RNAAlphabet{4}, 0b10000)
        end

        @testset "AminoAcid" begin
            @test decode(AminoAcidAlphabet, 0x00) === AA_A
            for x in 0x00:0x1b
                @test decode(AminoAcidAlphabet, x) === convert(AminoAcid, x)
            end
            @test_throws Seq.DecodeError decode(AminoAcidAlphabet, 0x1c)
        end
    end
end

@testset "BioSequence" begin
    @testset "Constructing empty sequences" begin
        @test DNASequence() == BioSequence(DNA)
        @test RNASequence() == BioSequence(RNA)
        @test AminoAcidSequence() == BioSequence(AminoAcid)
        @test CharSequence() == BioSequence(Char)
    end

    @testset "Constructing uninitialized sequences" begin
        @test isa(BioSequence{DNAAlphabet{2}}(0), BioSequence)
        @test isa(BioSequence{DNAAlphabet{4}}(10), BioSequence)
        @test isa(BioSequence{RNAAlphabet{2}}(0), BioSequence)
        @test isa(BioSequence{RNAAlphabet{4}}(10), BioSequence)
        @test isa(BioSequence{AminoAcidAlphabet}(10), BioSequence)
    end

    @testset "Conversion from/to strings" begin
        # Check that sequences in strings survive round trip conversion:
        #   String → BioSequence → String
        function test_string_construction(A::Type, seq::AbstractString)
            @test convert(AbstractString, BioSequence{A}(seq)) == uppercase(seq)
        end

        function test_string_parse(A::Type, seq::AbstractString)
            @test parse(BioSequence{A}, seq) == BioSequence{A}(seq)
        end

        for len in [0, 1, 2, 3, 10, 32, 1000, 10000]
            test_string_construction(DNAAlphabet{4}, random_dna(len))
            test_string_construction(DNAAlphabet{4}, lowercase(random_dna(len)))
            test_string_construction(RNAAlphabet{4}, lowercase(random_rna(len)))
            test_string_construction(RNAAlphabet{4}, random_rna(len))
            test_string_construction(AminoAcidAlphabet, random_aa(len))
            test_string_construction(AminoAcidAlphabet, lowercase(random_aa(len)))

            test_string_parse(DNAAlphabet{4}, random_dna(len))
            test_string_parse(DNAAlphabet{4}, lowercase(random_dna(len)))
            test_string_parse(RNAAlphabet{4}, lowercase(random_rna(len)))
            test_string_parse(RNAAlphabet{4}, random_rna(len))
            test_string_parse(AminoAcidAlphabet, random_aa(len))
            test_string_parse(AminoAcidAlphabet, lowercase(random_aa(len)))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_string_construction(DNAAlphabet{2}, random_dna(len, probs))
            test_string_construction(DNAAlphabet{2}, lowercase(random_dna(len, probs)))
            test_string_construction(RNAAlphabet{2}, random_rna(len, probs))
            test_string_construction(RNAAlphabet{2}, lowercase(random_rna(len, probs)))

            test_string_parse(DNAAlphabet{2}, random_dna(len, probs))
            test_string_parse(DNAAlphabet{2}, lowercase(random_dna(len, probs)))
            test_string_parse(RNAAlphabet{2}, random_rna(len, probs))
            test_string_parse(RNAAlphabet{2}, lowercase(random_rna(len, probs)))
        end

        # non-standard string literal
        @test isa(dna"ACGTMRWSYKVHDBN-", DNASequence)
        @test isa(rna"ACGUMRWSYKVHDBN-", RNASequence)
        @test isa(aa"ARNDCQEGHILKMFPSTWYVBJZXOU*-", AminoAcidSequence)
        @test isa(char"いろは αβγ 甲乙丙", CharSequence)

        # Non-nucleotide characters should throw
        @test_throws Exception DNASequence("ACCNNCATTTTTTAGATXATAG")
        @test_throws Exception RNASequence("ACCNNCATTTTTTAGATXATAG")
        @test_throws Exception AminoAcidSequence("ATGHLMY@ZACAGNM")
    end

    @testset "Construction from vectors" begin
        function test_vector_construction(A, seq::AbstractString)
            T = eltype(A)
            xs = T[convert(T, c) for c in seq]
            @test BioSequence{A}(xs) == BioSequence{A}(seq)
        end

        for len in [0, 1, 10, 32, 1000, 10000]
            test_vector_construction(DNAAlphabet{4}, random_dna(len))
            test_vector_construction(RNAAlphabet{4}, random_rna(len))
            test_vector_construction(AminoAcidAlphabet, random_aa(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_vector_construction(DNAAlphabet{2}, random_dna(len, probs))
            test_vector_construction(RNAAlphabet{2}, random_rna(len, probs))
        end
    end

    @testset "Conversion between 2-bit and 4-bit encodings" begin
        function test_conversion(A1, A2, seq)
            @test convert(BioSequence{A1}, BioSequence{A2}(seq)) == convert(BioSequence{A1}, seq)
        end

        test_conversion(DNAAlphabet{2}, DNAAlphabet{4}, "")
        test_conversion(DNAAlphabet{4}, DNAAlphabet{2}, "")
        test_conversion(RNAAlphabet{4}, RNAAlphabet{2}, "")
        test_conversion(RNAAlphabet{2}, RNAAlphabet{4}, "")

        test_conversion(DNAAlphabet{2}, DNAAlphabet{4}, "ACGT")
        test_conversion(DNAAlphabet{4}, DNAAlphabet{2}, "ACGT")
        test_conversion(RNAAlphabet{4}, RNAAlphabet{2}, "ACGU")
        test_conversion(RNAAlphabet{2}, RNAAlphabet{4}, "ACGU")

        test_conversion(DNAAlphabet{2}, DNAAlphabet{4}, "ACGT"^100)
        test_conversion(DNAAlphabet{4}, DNAAlphabet{2}, "ACGT"^100)
        test_conversion(RNAAlphabet{4}, RNAAlphabet{2}, "ACGU"^100)
        test_conversion(RNAAlphabet{2}, RNAAlphabet{4}, "ACGU"^100)

        # ambiguous nucleotides cannot be stored in 2-bit encoding
        EncodeError = Seq.EncodeError
        @test_throws EncodeError convert(BioSequence{DNAAlphabet{2}}, dna"AN")
        @test_throws EncodeError convert(BioSequence{RNAAlphabet{2}}, rna"AN")

        # test promotion
        a = BioSequence{DNAAlphabet{2}}("ATCG")
        b = BioSequence{DNAAlphabet{4}}("ATCG")
        c = BioSequence{RNAAlphabet{2}}("AUCG")
        d = BioSequence{RNAAlphabet{4}}("AUCG")

        @test typeof(promote(a, b)) == Tuple{BioSequence{DNAAlphabet{4}},BioSequence{DNAAlphabet{4}}}
        @test typeof(promote(c, d)) == Tuple{BioSequence{RNAAlphabet{4}},BioSequence{RNAAlphabet{4}}}
        @test typeof(promote(a, d)) == Tuple{BioSequence{DNAAlphabet{2}},BioSequence{RNAAlphabet{4}}}
        @test typeof(promote(a, b, d)) == Tuple{BioSequence{DNAAlphabet{2}},BioSequence{DNAAlphabet{4}},BioSequence{RNAAlphabet{4}}}
    end

    @testset "Conversion between RNA and DNA" begin
        @test convert(RNASequence, DNASequence("ACGTN")) == rna"ACGUN"
        @test convert(DNASequence, RNASequence("ACGUN")) == dna"ACGTN"
    end

    @testset "Conversion to Matrices" begin
        dna = [dna"AAA", dna"TTT", dna"CCC", dna"GGG"]
        dnathrow = [dna"AAA", dna"TTTAAA", dna"CCC", dna"GGG"]

        rna = [rna"AAA", rna"UUU", rna"CCC", rna"GGG"]
        rnathrow = [rna"AAA", rna"UUU", rna"CCCUUU", rna"GGG"]

        prot = [aa"AMG", aa"AMG", aa"AMG", aa"AMG"]

        sitemajdna = [
            DNA_A  DNA_A  DNA_A
            DNA_T  DNA_T  DNA_T
            DNA_C  DNA_C  DNA_C
            DNA_G  DNA_G  DNA_G
        ]
        seqmajdna = [
            DNA_A  DNA_T  DNA_C  DNA_G
            DNA_A  DNA_T  DNA_C  DNA_G
            DNA_A  DNA_T  DNA_C  DNA_G
        ]
        sitemajrna = [
            RNA_A  RNA_A  RNA_A
            RNA_U  RNA_U  RNA_U
            RNA_C  RNA_C  RNA_C
            RNA_G  RNA_G  RNA_G
        ]
        seqmajrna = [
            RNA_A  RNA_U  RNA_C  RNA_G
            RNA_A  RNA_U  RNA_C  RNA_G
            RNA_A  RNA_U  RNA_C  RNA_G
        ]
        sitemajnucint = [
            0x01  0x01  0x01
            0x08  0x08  0x08
            0x02  0x02  0x02
            0x04  0x04  0x04
        ]
        seqmajnucint = [
            0x01  0x08  0x02  0x04
            0x01  0x08  0x02  0x04
            0x01  0x08  0x02  0x04
        ]
        sitemajaa = [
            AA_A AA_M AA_G
            AA_A AA_M AA_G
            AA_A AA_M AA_G
            AA_A AA_M AA_G
        ]
        seqmajaa = [
            AA_A AA_A AA_A AA_A
            AA_M AA_M AA_M AA_M
            AA_G AA_G AA_G AA_G
        ]

        @test seqmatrix(dna, :site) == sitemajdna
        @test seqmatrix(rna, :site) == sitemajrna
        @test seqmatrix(prot, :site) == sitemajaa
        @test seqmatrix(UInt8, dna, :site) == sitemajnucint
        @test seqmatrix(UInt8, rna, :site) == sitemajnucint

        @test seqmatrix(dna, :seq) == seqmajdna
        @test seqmatrix(rna, :seq) == seqmajrna
        @test seqmatrix(prot, :seq) == seqmajaa
        @test seqmatrix(UInt8, dna, :seq) == seqmajnucint
        @test seqmatrix(UInt8, rna, :seq) == seqmajnucint

        @test seqmatrix([dna"", dna"", dna""], :site) == Matrix{DNA}(3, 0)
        @test seqmatrix([dna"", dna"", dna""], :seq) == Matrix{DNA}(0, 3)
        @test seqmatrix([rna"", rna"", rna""], :site) == Matrix{RNA}(3, 0)
        @test seqmatrix([rna"", rna"", rna""], :seq) == Matrix{RNA}(0, 3)
        @test seqmatrix(UInt8, [dna"", dna"", dna""], :site) == Matrix{UInt8}(3, 0)
        @test seqmatrix(UInt8, [dna"", dna"", dna""], :seq) == Matrix{UInt8}(0, 3)
        @test seqmatrix(UInt8, [rna"", rna"", rna""], :site) == Matrix{UInt8}(3, 0)
        @test seqmatrix(UInt8, [rna"", rna"", rna""], :seq) == Matrix{UInt8}(0, 3)

        @test_throws ArgumentError seqmatrix(dnathrow, :site)
        @test_throws ArgumentError seqmatrix(rnathrow, :seq)
        @test_throws ArgumentError seqmatrix(dna, :lol)
        @test_throws MethodError seqmatrix(AminoAcid, dna, :site)
        @test_throws ArgumentError seqmatrix(DNASequence[], :site)
        @test_throws ArgumentError seqmatrix(DNASequence[], :seq)

    end

    @testset "Copy" begin
        function test_copy(A, seq)
            @test convert(AbstractString, copy(BioSequence{A}(seq))) == seq
        end

        for len in [1, 10, 16, 32, 1000, 10000]
            test_copy(DNAAlphabet{4}, random_dna(len))
            test_copy(RNAAlphabet{4}, random_rna(len))
            test_copy(AminoAcidAlphabet, random_aa(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_copy(DNAAlphabet{2}, random_dna(len, probs))
            test_copy(RNAAlphabet{2}, random_rna(len, probs))
        end

        seq = dna"ACGTACGTACGT"
        subseq = seq[3:6]
        @test copy(subseq) == dna"GTAC"
    end

    @testset "Concatenation" begin
        function test_concatenation(A, chunks)
            parts = UnitRange{Int}[]
            for i in 1:endof(chunks)
                start = rand(1:length(chunks[i]))
                stop = rand(start:length(chunks[i]))
                push!(parts, start:stop)
            end
            str = string([chunk[parts[i]] for (i, chunk) in enumerate(chunks)]...)
            seq = *([BioSequence{A}(chunk)[parts[i]] for (i, chunk) in enumerate(chunks)]...)
            @test convert(AbstractString, seq) == uppercase(str)
        end

        for _ in 1:100
            n = rand(1:10)
            chunks = [random_dna(rand(1:100)) for _ in 1:n]
            test_concatenation(DNAAlphabet{4}, chunks)

            chunks = [random_rna(rand(1:100)) for _ in 1:n]
            test_concatenation(RNAAlphabet{4}, chunks)

            chunks = [random_aa(rand(1:100)) for _ in 1:n]
            test_concatenation(AminoAcidAlphabet, chunks)

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            chunks = [random_dna(rand(1:100), probs) for _ in 1:n]
            test_concatenation(DNAAlphabet{2}, chunks)

            chunks = [random_rna(rand(1:100), probs) for _ in 1:n]
            test_concatenation(RNAAlphabet{2}, chunks)
        end
    end

    @testset "Repetition" begin
        function test_repetition(A, chunk)
            start = rand(1:length(chunk))
            stop = rand(start:length(chunk))
            n = rand(1:10)
            str = chunk[start:stop] ^ n
            seq = BioSequence{A}(chunk)[start:stop] ^ n
            @test convert(AbstractString, seq) == uppercase(str)
        end

        for _ in 1:10
            chunk = random_dna(rand(1:100))
            test_repetition(DNAAlphabet{4}, chunk)

            chunk = random_rna(rand(1:100))
            test_repetition(RNAAlphabet{4}, chunk)

            chunk = random_aa(rand(1:100))
            test_repetition(AminoAcidAlphabet, chunk)

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            chunk = random_dna(rand(1:100), probs)
            test_repetition(DNAAlphabet{2}, chunk)

            chunk = random_rna(rand(1:100), probs)
            test_repetition(RNAAlphabet{2}, chunk)
        end
    end

    @testset "Equality" begin
        @testset "DNA" begin
            a = dna"ACTGN"
            b = dna"ACTGN"
            @test a == b
            @test dna"ACTGN" == dna"ACTGN"
            @test dna"ACTGN" != dna"ACTGA"
            @test dna"ACTGN" != dna"ACTG"
            @test dna"ACTG"  != dna"ACTGN"

            a = dna"ACGTNACGTN"
            b = dna"""
            ACGTN
            ACGTN
            """
            @test a == b
        end

        @testset "RNA" begin
            a = rna"ACUGN"
            b = rna"ACUGN"
            @test a == b
            @test rna"ACUGN" == rna"ACUGN"
            @test rna"ACUGN" != rna"ACUGA"
            @test rna"ACUGN" != rna"ACUG"
            @test rna"ACUG"  != rna"ACUGN"

            a = rna"ACUGNACUGN"
            b = rna"""
            ACUGN
            ACUGN
            """
            @test a == b
        end

        @testset "AminoAcid" begin
            a = aa"ARNDCQEGHILKMFPSTWYVX"
            b = aa"ARNDCQEGHILKMFPSTWYVX"
            @test a == b
            @test a != aa"ARNDCQEGHILKMFPSTWYXV"
            @test a != aa"ARNDCQEGHLKMFPSTWYVX"

            b = aa"""
            ARNDCQEGHI
            LKMFPSTWYV
            X
            """
            @test a == b
        end
    end

    @testset "Hash" begin
        seq = dna"ACGTACGT"
        @test isa(hash(seq), UInt64)
        @test hash(seq) === hash(dna"ACGTACGT")
        @test hash(seq) !== hash(seq[1:6])
        @test hash(seq) !== hash(seq[1:7])
        @test hash(seq) === hash(seq[1:8])
        @test hash(seq[1:4]) === hash(dna"ACGT")
        @test hash(seq[2:4]) === hash(dna"CGT")
        @test hash(seq[3:4]) === hash(dna"GT")
        @test hash(seq[4:4]) === hash(dna"T")
        @test hash(seq[5:8]) === hash(dna"ACGT")

        @test hash(dna"") !== hash(dna"A")
        @test hash(dna"A") !== hash(dna"AA")
        @test hash(dna"AA") !== hash(dna"AAA")
        @test hash(dna"AAA") !== hash(dna"AAAA")

        for n in 1:200, seq in [dna"A", dna"AC", dna"ACG", dna"ACGT", dna"ACGTN"]
            @test hash(seq^n) === hash((dna""     * seq^n)[1:end])
            @test hash(seq^n) === hash((dna"T"    * seq^n)[2:end])
            @test hash(seq^n) === hash((dna"TT"   * seq^n)[3:end])
            @test hash(seq^n) === hash((dna"TTT"  * seq^n)[4:end])
            @test hash(seq^n) === hash((dna"TTTT" * seq^n)[5:end])

            @test hash(seq^n) === hash((seq^n * dna""    )[1:end  ])
            @test hash(seq^n) === hash((seq^n * dna"T"   )[1:end-1])
            @test hash(seq^n) === hash((seq^n * dna"TT"  )[1:end-2])
            @test hash(seq^n) === hash((seq^n * dna"TTT" )[1:end-3])
            @test hash(seq^n) === hash((seq^n * dna"TTTT")[1:end-4])
        end

        @test hash(rna"AAUU") === hash(rna"AAUU")
        @test hash(rna"AAUUAA"[3:5]) === hash(rna"UUA")
        @test hash(aa"MTTQAPMFTQPLQ") === hash(aa"MTTQAPMFTQPLQ")
        @test hash(aa"MTTQAPMFTQPLQ"[5:10]) === hash(aa"APMFTQ")

        @testset "MinHash" begin
            seq = DNASequence(random_dna(1000))
            h = minhash(seq, 10, 100)

            @test length(h) == 100
            @test h == minhash(seq, 10, 100)

            @test_throws BoundsError h[101]
        end

    end

    @testset "Length" begin
        for len in [0, 1, 2, 3, 10, 16, 32, 1000, 10000]
            seq = DNASequence(random_dna(len))
            @test length(seq) === endof(seq) === len

            seq = RNASequence(random_rna(len))
            @test length(seq) === endof(seq) === len

            seq = AminoAcidSequence(random_aa(len))
            @test length(seq) === endof(seq) === len
        end

        @test length(char"いろはabc") === 6
    end

    @testset "Access" begin
        dna_seq = dna"ACTG"

        @test dna_seq[1] === DNA_A
        @test dna_seq[2] === DNA_C
        @test dna_seq[3] === DNA_T
        @test dna_seq[4] === DNA_G

        # Access indexes out of bounds
        @test_throws BoundsError dna_seq[-1]
        @test_throws BoundsError dna_seq[0]
        @test_throws BoundsError dna_seq[5]

        @test dna"ACTGNACTGN"[1:5] == dna"ACTGN"
        @test dna"ACTGNACTGN"[5:1] == dna""

        rna_seq = rna"ACUG"
        @test rna_seq[1] === RNA_A
        @test rna_seq[2] === RNA_C
        @test rna_seq[3] === RNA_U
        @test rna_seq[4] === RNA_G

        # Access indexes out of bounds
        @test_throws BoundsError rna_seq[-1]
        @test_throws BoundsError rna_seq[0]
        @test_throws BoundsError rna_seq[5]

        @test rna"ACUGNACUGN"[1:5] == rna"ACUGN"
        @test rna"ACUGNACUGN"[5:1] == rna""

        @test aa"KSAAV"[3] == AA_A
        @test char"いろはにほ"[3] == 'は'
    end

    @testset "Iteration" begin
        dna_seq = dna"ACTG"
        dna_vec = [DNA_A, DNA_C, DNA_T, DNA_G]
        @test all([nt === dna_vec[i] for (i, nt) in enumerate(dna_seq)])

        rna_seq = rna"ACUG"
        rna_vec = [RNA_A, RNA_C, RNA_U, RNA_G]
        @test all([nt === rna_vec[i] for (i, nt) in enumerate(rna_seq)])

        aa_seq = aa"ARNPS"
        aa_vec = [AA_A, AA_R, AA_N, AA_P, AA_S]
        @test all([aa == aa_vec[i] for (i, aa) in enumerate(aa_seq)])
    end

    @testset "Subsequence construction" begin
        function test_subseq(A, seq)
            bioseq = BioSequence{A}(seq)
            for _ in 1:100
                part = random_interval(1, endof(seq))
                @test convert(AbstractString, bioseq[part]) == seq[part]
            end
        end

        for len in [1, 10, 32, 1000, 10000, 100000]
            test_subseq(DNAAlphabet{4}, random_dna(len))
            test_subseq(RNAAlphabet{4}, random_rna(len))
            test_subseq(AminoAcidAlphabet, random_aa(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_subseq(DNAAlphabet{2}, random_dna(len, probs))
            test_subseq(RNAAlphabet{2}, random_rna(len, probs))
        end

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

    @testset "Mutability" begin
        @testset "setindex!" begin
            s = dna"ACGT"
            s[1] = DNA_A
            @test s == dna"ACGT"
            s[1] = DNA_C
            @test s == dna"CCGT"
            s[2] = DNA_T
            @test s == dna"CTGT"
            s[4] = DNA_A
            @test s == dna"CTGA"
            @test_throws BoundsError s[0]
            @test_throws BoundsError s[5]

            s = dna"ACGTACGT"
            s[3:5] = DNA_A
            @test s == dna"ACAAACGT"
            s[3:5] = DNA_T
            @test s == dna"ACTTTCGT"
            s[1:8] = DNA_C
            @test s == dna"CCCCCCCC"
            s[:] = DNA_G
            @test s == dna"GGGGGGGG"
            @test_throws BoundsError s[0:3]
            @test_throws BoundsError s[5:10]

            s = dna"ACGTACGT"
            s[[1,2,3]] = DNA_T
            @test s == dna"TTTTACGT"
            s[[1,5,8]] = DNA_C
            @test s == dna"CTTTCCGC"
            @test_throws BoundsError s[[3,9]] = DNA_A

            s = dna"ACGT"
            s[[true, false, false, true]] = DNA_G
            @test s == dna"GCGG"
            s[trues(4)] = DNA_A
            @test s == dna"AAAA"
            @test_throws BoundsError s[[true, false, false]] = DNA_G
            @test_throws BoundsError s[[true, false, false, false, true]] = DNA_G

            s = dna"ACGTACGT"
            s[2:3] = dna"AA"
            @test s == dna"AAATACGT"
            s[7:8] = dna"CC"
            @test s == dna"AAATACCC"
            s[:] = dna"AACCGGTT"
            @test s == dna"AACCGGTT"
            @test_throws BoundsError       s[0:1] = dna"AA"
            @test_throws DimensionMismatch s[3:4] = dna"A"

            s = dna"ACGTACGT"
            s[[1,4]] = dna"TA"
            @test s == dna"TCGAACGT"
            s[[2,3,5]] = dna"CAT"
            @test s == dna"TCAATCGT"
            @test_throws BoundsError       s[[1,2,9]] = dna"AAA"
            @test_throws DimensionMismatch s[[1,2,8]] = dna"AA"

            s = dna"ACGT"
            s[[true,false,true,false]] = dna"TT"
            @test s == dna"TCTT"
            s[trues(4)] = dna"AAAA"
            @test s == dna"AAAA"
            @test_throws BoundsError       s[[true,false,true]] = dna"TT"
            @test_throws DimensionMismatch s[[true,false,true,true]] = dna"TT"
        end

        @testset "resize!" begin
            seq = dna""
            resize!(seq, 100)
            @test length(seq) == 100
            resize!(seq, 200)
            @test length(seq) == 200
            resize!(seq,  10)
            @test length(seq) == 10
            @test_throws ArgumentError resize!(seq, -1)
        end

        @testset "empty!" begin
            seq = dna"ACG"
            @test empty!(seq) == dna""
            @test length(seq) == 0
        end

        @testset "push!" begin
            seq = dna""
            @test push!(seq, DNA_A) == dna"A"
            @test push!(seq, DNA_C) == dna"AC"
            @test seq == dna"AC"
        end

        @testset "unshift!" begin
            seq = dna""
            @test unshift!(seq, DNA_A) == dna"A"
            @test unshift!(seq, DNA_C) == dna"CA"
            @test seq == dna"CA"
        end

        @testset "pop!" begin
            seq = dna"ACGT"
            @test pop!(seq) === DNA_T
            @test seq == dna"ACG"
            @test pop!(seq) === DNA_G
            @test seq == dna"AC"
            @test_throws ArgumentError pop!(dna"")
        end

        @testset "shift!" begin
            seq = dna"ACGT"
            @test shift!(seq) === DNA_A
            @test seq == dna"CGT"
            @test shift!(seq) === DNA_C
            @test seq == dna"GT"
            @test_throws ArgumentError shift!(dna"")
        end

        @testset "insert!" begin
            seq = dna"ACGT"
            @test insert!(seq, 2, DNA_G) == dna"AGCGT"
            @test insert!(seq, 5, DNA_A) == dna"AGCGAT"
            @test_throws BoundsError insert!(seq, 10, DNA_T)
        end

        @testset "deleteat!" begin
            seq = dna"ACGT"
            @test deleteat!(seq, 1) == dna"CGT"
            @test deleteat!(seq, 2) == dna"CT"
            @test_throws BoundsError deleteat!(seq, 10)

            seq = dna"ACGTACGT"
            @test deleteat!(seq, 3:5) == dna"ACCGT"
            @test_throws BoundsError deleteat!(seq, 10:12)
        end

        @testset "append!" begin
            seq = dna""
            @test append!(seq, dna"A") == dna"A"
            @test append!(seq, dna"ACG") == dna"AACG"
        end

        @testset "copy!" begin
            seq = dna"GGG"
            @test copy!(seq, dna"ACG") == dna"ACG"
            @test copy!(seq, dna"TTA") == dna"TTA"

            seq = dna"TCCC"
            @test copy!(seq, 2, dna"TT", 1) == dna"TTTC"
            seq = dna"TCCC"
            @test copy!(seq, 2, dna"TT", 1, 1) == dna"TTCC"

            seq = dna"ACGT"
            @test copy!(seq, seq) == dna"ACGT"
            @test copy!(seq, 1, seq, 3, 2) == dna"GTGT"
            seq = dna"ACGT"
            @test copy!(seq, 3, seq, 1, 2) == dna"ACAC"
        end

        @testset "orphan!" begin
            seq = repeat(dna"ACGT", 8)
            subseq = seq[16:17]
            Seq.orphan!(subseq)
            @test subseq == dna"TA"
        end
    end

    @testset "Print" begin
        buf = IOBuffer()
        print(buf, dna"")
        @test takebuf_string(buf) == ""

        buf = IOBuffer()
        print(buf, dna"ACGTN")
        @test takebuf_string(buf) == "ACGTN"

        buf = IOBuffer()
        print(buf, rna"ACGUN")
        @test takebuf_string(buf) == "ACGUN"

        buf = IOBuffer()
        print(buf, dna"A"^100)
        @test takebuf_string(buf) == "A"^100

        buf = IOBuffer()
        print(buf, dna"A"^100, width=70)
        @test takebuf_string(buf) == string("A"^70, '\n', "A"^30)
    end

    @testset "Transformations" begin
        function test_reverse(A, seq)
            revseq = reverse(BioSequence{A}(seq))
            @test convert(AbstractString, revseq) == reverse(seq)
        end

        function test_dna_complement(A, seq)
            comp = complement(BioSequence{A}(seq))
            @test convert(AbstractString, comp) == dna_complement(seq)
        end

        function test_rna_complement(A, seq)
            comp = complement(BioSequence{A}(seq))
            @test convert(AbstractString, comp) == rna_complement(seq)
        end

        function test_dna_revcomp(A, seq)
            revcomp = reverse_complement(BioSequence{A}(seq))
            @test convert(AbstractString, revcomp) == reverse(dna_complement(seq))
        end

        function test_rna_revcomp(A, seq)
            revcomp = reverse_complement(BioSequence{A}(seq))
            @test convert(AbstractString, revcomp) == reverse(rna_complement(seq))
        end

        @testset "Reverse" begin
            for len in [0, 1, 10, 32, 1000, 10000, 100000], _ in 1:10
                test_reverse(DNAAlphabet{4}, random_dna(len))
                test_reverse(RNAAlphabet{4}, random_rna(len))
                test_reverse(AminoAcidAlphabet, random_aa(len))

                probs = [0.25, 0.25, 0.25, 0.25, 0.00]
                test_reverse(DNAAlphabet{2}, random_dna(len, probs))
                test_reverse(RNAAlphabet{2}, random_rna(len, probs))
            end
        end

        @testset "Complement" begin
            for len in [0, 1, 10, 32, 1000, 10000, 100000], _ in 1:10
                test_dna_complement(DNAAlphabet{4}, random_dna(len))
                test_rna_complement(RNAAlphabet{4}, random_rna(len))

                probs = [0.25, 0.25, 0.25, 0.25, 0.00]
                test_dna_complement(DNAAlphabet{2}, random_dna(len, probs))
                test_rna_complement(RNAAlphabet{2}, random_rna(len, probs))
            end
        end

        @testset "Reverse complement" begin
            for len in [0, 1, 10, 32, 1000, 10000, 100000], _ in 1:10
                test_dna_revcomp(DNAAlphabet{4}, random_dna(len))
                test_rna_revcomp(RNAAlphabet{4}, random_rna(len))

                probs = [0.25, 0.25, 0.25, 0.25, 0.00]
                test_dna_revcomp(DNAAlphabet{2}, random_dna(len, probs))
                test_rna_revcomp(RNAAlphabet{2}, random_rna(len, probs))
            end
        end

        @testset "Map" begin
            seq = dna""
            @test map(identity, seq) == dna""
            seq = dna"AAA"
            @test map(x -> DNA_C, seq) == dna"CCC"
            seq = dna"ACGTNACGTN"
            @test map(x -> x == DNA_N ? DNA_A : x, seq) == dna"ACGTAACGTA"
            @test seq == dna"ACGTNACGTN"
            @test map!(x -> x == DNA_N ? DNA_A : x, seq) === seq
            @test seq == dna"ACGTAACGTA"

            seq = rna""
            @test map(identity, seq) == rna""
            seq = rna"AAA"
            @test map(x -> RNA_C, seq) == rna"CCC"
            seq = rna"ACGUNACGUN"
            @test map(x -> x == RNA_N ? RNA_A : x, seq) == rna"ACGUAACGUA"
            @test seq == rna"ACGUNACGUN"
            @test map!(x -> x == RNA_N ? RNA_A : x, seq) === seq
            @test seq == rna"ACGUAACGUA"

            seq = aa""
            @test map(identity, seq) == aa""
            seq = aa"MMM"
            @test map(x -> AA_P, seq) == aa"PPP"
            seq = aa"XRNDCQXE"
            @test map(x -> x == AA_X ? AA_A : x, seq) == aa"ARNDCQAE"
            @test seq == aa"XRNDCQXE"
            @test map!(x -> x == AA_X ? AA_A : x, seq) === seq
            @test seq == aa"ARNDCQAE"
        end

        @testset "Filter" begin
            seq = dna""
            @test filter(x -> true, seq) == dna""
            @test filter(x -> false, seq) == dna""
            seq = dna"AAA"
            @test filter(x -> x == DNA_A, seq) == dna"AAA"
            @test filter(x -> x == DNA_C, seq) == dna""
            seq = dna"ACGTNACGTN"
            @test filter(x -> x == DNA_N, seq) == dna"NN"
            @test filter(x -> x != DNA_N, seq) == dna"ACGTACGT"
            @test seq == dna"ACGTNACGTN"
            @test filter!(x -> x != DNA_N, seq) == seq
            @test seq == dna"ACGTACGT"

            for len in 1:50, _ in 1:10
                str = random_dna(len)
                seq = DNASequence(str)
                @test filter(x -> x == DNA_N, seq) ==
                    DNASequence(filter(x -> x == 'N', str))
                @test filter(x -> x != DNA_N, seq) ==
                    DNASequence(filter(x -> x != 'N', str))
            end

            seq = rna""
            @test filter(x -> true, seq) == rna""
            @test filter(x -> false, seq) == rna""
            seq = rna"AAA"
            @test filter(x -> x == RNA_A, seq) == rna"AAA"
            @test filter(x -> x == RNA_C, seq) == rna""
            seq = rna"ACGUNACGUN"
            @test filter(x -> x == RNA_N, seq) == rna"NN"
            @test filter(x -> x != RNA_N, seq) == rna"ACGUACGU"
            @test seq == rna"ACGUNACGUN"
            @test filter!(x -> x != RNA_N, seq) == seq
            @test seq == rna"ACGUACGU"

            seq = aa""
            @test filter(x -> true, seq) == aa""
            @test filter(x -> false, seq) == aa""
            seq = aa"PPP"
            @test filter(x -> x == AA_P, seq) == aa"PPP"
            @test filter(x -> x != AA_P, seq) == aa""
            seq = aa"ARNDCQXGHILKMFPXTWYVOUX"
            @test filter(x -> x == AA_X, seq) == aa"XXX"
            @test filter(x -> x != AA_X, seq) == aa"ARNDCQGHILKMFPTWYVOU"
            @test seq == aa"ARNDCQXGHILKMFPXTWYVOUX"
            @test filter!(x -> x != AA_X, seq) == seq
            @test seq == aa"ARNDCQGHILKMFPTWYVOU"
        end
    end

    @testset "Predicates" begin
        # ispalindromic
        @test  ispalindromic(dna"")
        @test !ispalindromic(dna"A")
        @test !ispalindromic(dna"C")
        @test  ispalindromic(dna"AT")
        @test  ispalindromic(dna"CG")
        @test !ispalindromic(dna"AC")
        @test !ispalindromic(dna"TT")
        @test  ispalindromic(dna"ANT")
        @test  ispalindromic(dna"ACGT")
        @test !ispalindromic(dna"ACNT")

        @test  ispalindromic(DNAKmer("ACGT"))
        @test !ispalindromic(DNAKmer("CACG"))

        @test  ispalindromic(rna"")
        @test !ispalindromic(rna"A")
        @test !ispalindromic(rna"C")
        @test  ispalindromic(rna"AU")
        @test  ispalindromic(rna"CG")
        @test !ispalindromic(rna"AC")
        @test !ispalindromic(rna"UU")
        @test  ispalindromic(rna"ANU")
        @test  ispalindromic(rna"ACGU")
        @test !ispalindromic(rna"ACNU")

        @test  ispalindromic(RNAKmer("ACGU"))
        @test !ispalindromic(RNAKmer("CACG"))

        @test_throws Exception ispalindromic(aa"PQ")

        # hasambiguity
        @test !hasambiguity(dna"")
        @test !hasambiguity(dna"A")
        @test  hasambiguity(dna"N")
        @test !hasambiguity(dna"ACGT")
        @test  hasambiguity(dna"ANGT")

        @test !hasambiguity(DNAKmer("ACGT"))

        @test !hasambiguity(rna"")
        @test !hasambiguity(rna"A")
        @test  hasambiguity(rna"N")
        @test !hasambiguity(rna"ACGU")
        @test  hasambiguity(rna"ANGU")

        @test !hasambiguity(RNAKmer("ACGU"))

        @test !hasambiguity(aa"")
        @test !hasambiguity(aa"A")
        @test !hasambiguity(aa"P")
        @test  hasambiguity(aa"B")
        @test  hasambiguity(aa"X")
        @test !hasambiguity(aa"ARNDCQEGHILKMFPSTWYVOU")
        @test  hasambiguity(aa"ARXDCQEGHILKMFPSTWYVOU")

        # isrepetitive
        @test  isrepetitive(dna"")
        @test !isrepetitive(dna"", 1)
        @test  isrepetitive(dna"A")
        @test  isrepetitive(dna"A", 1)
        @test  isrepetitive(dna"AAA")
        @test !isrepetitive(dna"ACGT", 2)
        @test  isrepetitive(dna"AAGT", 2)
        @test  isrepetitive(dna"ACCG", 2)
        @test  isrepetitive(dna"ACGG", 2)
        @test !isrepetitive(dna"ACGTCCGT", 3)
        @test  isrepetitive(dna"ACGCCCGT", 3)

        @test !isrepetitive(DNAKmer("ACGT"), 2)
        @test  isrepetitive(DNAKmer("ACCT"), 2)
        @test !isrepetitive(DNAKmer("ACCT"), 3)
        @test  isrepetitive(DNAKmer("ACCCT"), 2)
        @test  isrepetitive(DNAKmer("ACCCT"), 3)

        @test  isrepetitive(rna"")
        @test !isrepetitive(rna"", 1)
        @test  isrepetitive(rna"A")
        @test  isrepetitive(rna"A", 1)
        @test  isrepetitive(rna"AAA")
        @test !isrepetitive(rna"ACGU", 2)
        @test  isrepetitive(rna"AAGU", 2)
        @test  isrepetitive(rna"ACCG", 2)
        @test  isrepetitive(rna"ACGG", 2)
        @test !isrepetitive(rna"ACGUCCGU", 3)
        @test  isrepetitive(rna"ACGCCCGU", 3)

        @test !isrepetitive(RNAKmer("ACGU"), 2)
        @test  isrepetitive(RNAKmer("ACCU"), 2)
        @test !isrepetitive(RNAKmer("ACCU"), 3)
        @test  isrepetitive(RNAKmer("ACCCU"), 2)
        @test  isrepetitive(RNAKmer("ACCCU"), 3)

        @test  isrepetitive(aa"")
        @test !isrepetitive(aa"PGQQ")
        @test  isrepetitive(aa"PGQQ", 2)
        @test !isrepetitive(aa"PPQQ", 3)
        @test  isrepetitive(aa"PPPQQ", 3)
    end

    @testset "Find" begin
        seq = dna"ACGNA"
        @test findnext(seq, DNA_A, 1) == 1
        @test findnext(seq, DNA_C, 1) == 2
        @test findnext(seq, DNA_G, 1) == 3
        @test findnext(seq, DNA_N, 1) == 4
        @test findnext(seq, DNA_T, 1) == 0
        @test findnext(seq, DNA_A, 2) == 5

        @test_throws BoundsError findnext(seq, DNA_A, 0)
        @test_throws BoundsError findnext(seq, DNA_A, 6)

        @test findprev(seq, DNA_A, 4) == 1
        @test findprev(seq, DNA_C, 4) == 2
        @test findprev(seq, DNA_G, 4) == 3
        @test findprev(seq, DNA_N, 4) == 4
        @test findprev(seq, DNA_T, 4) == 0
        @test findprev(seq, DNA_G, 2) == 0

        @test_throws BoundsError findprev(seq, DNA_A, 0)
        @test_throws BoundsError findprev(seq, DNA_A, 6)

        seq = dna"ACGNAN"
        @test findfirst(seq, DNA_A) == 1
        @test findfirst(seq, DNA_N) == 4
        @test findfirst(seq, DNA_T) == 0

        @test findlast(seq, DNA_A) == 5
        @test findlast(seq, DNA_N) == 6
        @test findlast(seq, DNA_T) == 0
    end

    @testset "Site counting" begin
        @testset "Naive methods" begin
            NC = Seq.NaiveCount

            alphabets = (DNAAlphabet{4}, DNAAlphabet{2},
                         RNAAlphabet{4}, RNAAlphabet{2})

            for alph in alphabets

                # Answers to these tests were worked out manually to verify
                # count_sites_naive was working correctly.
                # seqA and seqB contain all possible observations of sites.

                istwobit = Seq.bitsof(alph) == 2

                seqA, seqB = generate_possibilities_tester(alph)

                # Test methods which work on single sequences.
                @test count(Certain, NC, seqA) == ifelse(istwobit, length(seqA), 49)
                @test count(Certain, NC, seqB) == ifelse(istwobit, length(seqB), 19)
                @test count(Gap, NC, seqA) == ifelse(istwobit, 0, 16)
                @test count(Gap, NC, seqB) == ifelse(istwobit, 0, 1)
                @test count(Ambiguous, NC, seqA) == ifelse(istwobit, 0, length(seqA) - 65)
                @test count(Ambiguous, NC, seqB) == ifelse(istwobit, 0, length(seqB) - 20)

                # Test methods which work on two sequences.
                # Test when sequences are of the same bitencoding.

                @test count(Certain, NC, seqA, seqB) == count(Certain, NC, seqB, seqA) == 10
                @test count(Gap, NC, seqA, seqB) == count(Gap, NC, seqB, seqA) == ifelse(istwobit, 0, 16)
                @test count(Ambiguous, NC, seqA, seqB) == count(Ambiguous, NC, seqB, seqA) == ifelse(istwobit, 0, 121)
                @test count(Match, NC, seqA, seqB) == count(Match, NC, seqB, seqA) == length(alphabet(alph))
                @test count(Mismatch, NC, seqA, seqB) == count(Mismatch, NC, seqB, seqA) == (length(seqA) - length(alphabet(alph)))
            end

            # Test for when sequences are of different bitencodings.
            for alphs in [(DNAAlphabet{2}, DNAAlphabet{4}),
                          (RNAAlphabet{2}, RNAAlphabet{4})]
                seqA, seqB = generate_possibilities_tester(alphs...)
                @test count(Certain, NC, seqA, seqB) == count(Certain, NC, seqB, seqA) == 16
                @test count(Gap, NC, seqA, seqB) == count(Gap, NC, seqB, seqA) == 4
                @test count(Ambiguous, NC, seqA, seqB) == count(Ambiguous, NC, seqB, seqA) == 44
                @test count(Match, NC, seqA, seqB) == count(Match, NC, seqB, seqA) == 4
                @test count(Mismatch, NC, seqA, seqB) == count(Mismatch, NC, seqB, seqA) == 60
            end
        end

        @testset "Bit parallel methods" begin
            BC = Seq.BitparCount
            NC = Seq.NaiveCount
            # Having determined that naive counting algorithm is correct.
            # We can verify the bitparallel algorithm is correct by randomly
            # generating test-cases, and verifying the naive algorithm, and
            # the bitparallel algorithm give the same answer.
            @testset "4 bit encoding" begin
                alphabets = (DNAAlphabet{4}, RNAAlphabet{4})
                for alph in alphabets
                    for _ in 1:50
                        seqA = random_seq(alph, rand(10:100))
                        seqB = random_seq(alph, rand(10:100))
                        subA = seqA[1:rand(10:length(seqA))]
                        subB = seqB[1:rand(10:length(seqB))]
                        @test count(Mismatch, BC, subA, subB) == count(Mismatch, BC, subB, subA) == count(Mismatch, NC, subA, subB)
                        @test count(Match, BC, subA, subB) == count(Match, BC, subB, subA) == count(Match, NC, subA, subB)
                        @test count(Certain, BC, subA, subB) == count(Certain, BC, subB, subA) == count(Certain, NC, subA, subB)
                        @test count(Gap, BC, subA, subB) == count(Gap, BC, subB, subA) == count(Gap, NC, subA, subB)
                        @test count(Ambiguous, BC, subA, subB) == count(Ambiguous, BC, subB, subA) == count(Ambiguous, NC, subA, subB)
                    end
                end
            end

            @testset "2 bit encoding" begin
                alphabets = (DNAAlphabet{4}, RNAAlphabet{4})
                for alph in alphabets
                    for _ in 1:50
                        seqA = random_seq(alph, rand(10:100))
                        seqB = random_seq(alph, rand(10:100))
                        subA = seqA[1:rand(10:length(seqA))]
                        subB = seqB[1:rand(10:length(seqB))]
                        @test count(Mismatch, BC, subA, subB) == count(Mismatch, BC, subB, subA) == count(Mismatch, NC, subA, subB)
                        @test count(Match, BC, subA, subB) == count(Match, BC, subB, subA) == count(Match, NC, subA, subB)
                    end
                end
            end
        end

        @testset "Windowed methods" begin
            dnaA = dna"ATCGCCA-M"
            dnaB = dna"ATCGCCTAA"
            rnaA = rna"AUCGCCA-M"
            rnaB = rna"AUCGCCUAA"
            for seqs in ((dnaA, dnaB), (rnaA, rnaB))

                @test count(Certain, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                                 IntervalValue(2, 4, 3),
                                                                 IntervalValue(3, 5, 3),
                                                                 IntervalValue(4, 6, 3),
                                                                 IntervalValue(5, 7, 3),
                                                                 IntervalValue(6, 8, 2),
                                                                 IntervalValue(7, 9, 1)]
                @test count(Ambiguous, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                   IntervalValue(2, 4, 0),
                                                                   IntervalValue(3, 5, 0),
                                                                   IntervalValue(4, 6, 0),
                                                                   IntervalValue(5, 7, 0),
                                                                   IntervalValue(6, 8, 0),
                                                                   IntervalValue(7, 9, 1)]
                @test count(Gap, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                             IntervalValue(2, 4, 0),
                                                             IntervalValue(3, 5, 0),
                                                             IntervalValue(4, 6, 0),
                                                             IntervalValue(5, 7, 0),
                                                             IntervalValue(6, 8, 1),
                                                             IntervalValue(7, 9, 1)]
                @test count(Match, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                               IntervalValue(2, 4, 3),
                                                               IntervalValue(3, 5, 3),
                                                               IntervalValue(4, 6, 3),
                                                               IntervalValue(5, 7, 2),
                                                               IntervalValue(6, 8, 1),
                                                               IntervalValue(7, 9, 0)]
                @test count(Mismatch, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                  IntervalValue(2, 4, 0),
                                                                  IntervalValue(3, 5, 0),
                                                                  IntervalValue(4, 6, 0),
                                                                  IntervalValue(5, 7, 1),
                                                                  IntervalValue(6, 8, 2),
                                                                  IntervalValue(7, 9, 3)]
            end
        end

        @testset "Pairwise methods" begin
            dnas = [dna"ATCGCCA-", dna"ATCGCCTA", dna"ATCGCCT-", dna"GTCGCCTA"]
            rnas = [rna"AUCGCCA-", rna"AUCGCCUA", rna"AUCGCCU-", rna"GUCGCCUA"]
            answer_mismatch = PWM{Int, false}([0 2 1 3; 2 0 1 1; 1 1 0 2; 3 1 2 0])
            answer_match = PWM{Int, false}([0 6 7 5; 6 0 7 7; 7 7 0 6; 5 7 6 0])
            for i in (dnas, rnas)
                @test count_pairwise(Mismatch, i...) == answer_mismatch
                @test count_pairwise(Match, i...) == answer_match
                @test count_pairwise(Certain, i...) == PWM{Int, false}([0 7 7 7; 7 0 7 8; 7 7 0 7; 7 8 7 0])
                @test count_pairwise(Ambiguous, i...) == PWM{Int, false}([0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0])
                @test count_pairwise(Gap, i...) == PWM{Int, false}([0 1 1 1; 1 0 1 0; 1 1 0 1; 1 0 1 0])
            end
        end
    end

    @testset "GC content" begin
        @test gc_content(dna"") === 0.0
        @test gc_content(dna"AATA") === 0.0
        @test gc_content(dna"ACGT") === 0.5
        @test gc_content(dna"CGGC") === 1.0
        @test gc_content(dna"ACATTGTGTATAACAAAAGG") === 6 / 20

        @test gc_content(DNAKmer("")) === 0.0
        @test gc_content(DNAKmer("AATA")) === 0.0
        @test gc_content(DNAKmer("ACGT")) === 0.5
        @test gc_content(DNAKmer("CGGC")) === 1.0
        @test gc_content(DNAKmer("ACATTGTGTATAACAAAAGG")) === 6 / 20

        @test gc_content(rna"") === 0.0
        @test gc_content(rna"AAUA") === 0.0
        @test gc_content(rna"ACGU") === 0.5
        @test gc_content(rna"CGGC") === 1.0
        @test gc_content(rna"ACAUUGUGUAUAACAAAAGG") === 6 / 20

        @test gc_content(RNAKmer("")) === 0.0
        @test gc_content(RNAKmer("AAUA")) === 0.0
        @test gc_content(RNAKmer("ACGU")) === 0.5
        @test gc_content(RNAKmer("CGGC")) === 1.0
        @test gc_content(RNAKmer("ACAUUGUGUAUAACAAAAGG")) === 6 / 20

        @test_throws Exception gc_content(aa"ARN")
    end

    @testset "Ambiguous nucleotide positions" begin
        function test_ns(A, seq)
            expected = Int[]
            for i in 1:length(seq)
                if seq[i] == 'N'
                    push!(expected, i)
                end
            end
            bioseq = BioSequence{A}(seq)
            @test collect(ambiguous_positions(bioseq)) == expected
        end

        for len in [1, 10, 32, 1000, 10000, 100000]
            test_ns(DNAAlphabet{4}, random_dna(len))
            test_ns(RNAAlphabet{4}, random_rna(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_ns(DNAAlphabet{2}, random_dna(len, probs))
            test_ns(RNAAlphabet{2}, random_rna(len, probs))
        end

        dna_seq = dna"ANANANA"
        pos = [2, 4, 6]
        for (i, p) in enumerate(ambiguous_positions(dna_seq))
            @test p == pos[i]
        end

        dna_seq = dna"NATTCGRATY"
        pos = [1, 7, 10]
        for (i, p) in enumerate(ambiguous_positions(dna_seq))
            @test p == pos[i]
        end

        rna_seq = rna"ANANANA"
        pos = [2, 4, 6]
        for (i, p) in enumerate(ambiguous_positions(rna_seq))
            @test p == pos[i]
        end

        rna_seq = rna"NAUUCGRAUY"
        pos = [1, 7, 10]
        for (i, p) in enumerate(ambiguous_positions(rna_seq))
            @test p == pos[i]
        end
    end

    @testset "Shuffle" begin
        @testset for _ in 1:10
            @test shuffle(dna"") == dna""
            @test shuffle(dna"A") == dna"A"
            @test shuffle(dna"C") == dna"C"
        end

        seq = dna"ACGTN"^10
        @test shuffle(seq) != dna"ACGTN"^10
        @test seq == dna"ACGTN"^10
        @test shuffle!(seq) === seq
        @test seq != dna"ACGTN"^10
        @test count(x -> x == DNA_A, seq) == 10
        @test count(x -> x == DNA_C, seq) == 10
        @test count(x -> x == DNA_G, seq) == 10
        @test count(x -> x == DNA_T, seq) == 10
        @test count(x -> x == DNA_N, seq) == 10
    end
end

@testset "ReferenceSequence" begin
    @testset "Construction" begin
        @test isa(ReferenceSequence(dna""), ReferenceSequence)
        @test isa(ReferenceSequence(dna"ACGTN"^30), ReferenceSequence)
        @test isa(ReferenceSequence("ACGTN"^30), ReferenceSequence)
        @test_throws Exception ReferenceSequence(dna"ACGRT")
        @test_throws Exception ReferenceSequence("ACGRT")
    end

    @testset "Conversion" begin
        seq = DNASequence(random_dna(100))
        refseq = ReferenceSequence(seq)
        @test refseq == seq
        @test convert(DNASequence, refseq) == seq
    end

    @testset "Basic Operations" begin
        seq = ReferenceSequence(dna"")
        @test length(seq) === endof(seq) === 0
        @test isempty(seq)
        @test ReferenceSequence(dna"") == seq
        @test_throws BoundsError seq[0]
        @test_throws BoundsError seq[1]
        @test_throws BoundsError seq[2:3]
        @test collect(seq) == DNA[]

        seq = ReferenceSequence(dna"ACGTN")
        @test length(seq) === endof(seq) === 5
        @test !isempty(seq)
        @test seq[1] === DNA_A
        @test seq[2] === DNA_C
        @test seq[3] === DNA_G
        @test seq[4] === DNA_T
        @test seq[5] === DNA_N
        @test seq[2:3] == dna"CG"
        @test seq[3:5] == dna"GTN"
        @test seq[5:4] == dna""
        @test_throws BoundsError seq[0]
        @test_throws BoundsError seq[6]
        @test_throws BoundsError seq[0:2]
        @test_throws BoundsError seq[5:6]
        @test collect(seq) == [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
    end

    @testset "Long Sequence" begin
        str = random_dna(50000)
        @test ReferenceSequence(str) == DNASequence(str)
        str = ("N"^300 * random_dna(1000))^5 * "N"^300
        @test ReferenceSequence(str) == DNASequence(str)
    end

    @testset "Print" begin
        buf = IOBuffer()
        print(buf, ReferenceSequence(dna""))
        @test takebuf_string(buf) == ""

        buf = IOBuffer()
        print(buf, ReferenceSequence(dna"ACGTN"))
        @test takebuf_string(buf) == "ACGTN"

        buf = IOBuffer()
        print(buf, ReferenceSequence(dna"A"^100))
        @test takebuf_string(buf) == "A"^100
    end

    @testset "Random sequence" begin
        for _ in 1:10
            @test randdnaseq(0) == dna""
            @test randrnaseq(0) == rna""
            @test randaaseq(0) == aa""
        end

        for _ in 1:100
            len = rand(0:1000)
            @test length(randdnaseq(len)) == len
            @test length(randrnaseq(len)) == len
            @test length(randaaseq(len)) == len
        end

        for _ in 1:100
            len = 10_000

            counts = countmap(collect(randdnaseq(len)))
            counts_sum = 0
            for nt in dna"ACGT"
                # cdf(Binomial(10_000, 0.25), 2300) ≈ 1.68e-6
                @test 2300 < counts[nt]
                counts_sum += counts[nt]
            end
            @test counts_sum == len

            counts = countmap(collect(randrnaseq(len)))
            counts_sum = 0
            for nt in rna"ACGU"
                # cdf(Binomial(10_000, 0.25), 2300) ≈ 1.68e-6
                @test 2300 < counts[nt]
                counts_sum += counts[nt]
            end
            @test counts_sum == len

            counts = countmap(collect(randaaseq(len)))
            counts_sum = 0
            for aa in aa"ARNDCQEGHILKMFPSTWYV"
                # cdf(Binomial(10000, 0.05), 400) ≈ 1.21e-6
                @test 400 < counts[aa]
                counts_sum += counts[aa]
            end
            @test counts_sum == len
        end
    end
end

@testset "Composition" begin
    function string_nucleotide_count(::Type{DNA}, seq::AbstractString)
        counts = Dict{DNA, Int}(
            DNA_A => 0,
            DNA_C => 0,
            DNA_G => 0,
            DNA_T => 0,
            DNA_N => 0 )
        for c in seq
            counts[convert(DNA, c)] += 1
        end
        return counts
    end

    function string_nucleotide_count(::Type{RNA}, seq::AbstractString)
        counts = Dict{RNA, Int}(
            RNA_A => 0,
            RNA_C => 0,
            RNA_G => 0,
            RNA_U => 0,
            RNA_N => 0 )
        for c in seq
            counts[convert(RNA, c)] += 1
        end
        return counts
    end

    function check_nucleotide_count(::Type{DNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(DNA, seq)
        seq_counts = composition(DNASequence(seq))
        return string_counts[DNA_A] == seq_counts[DNA_A] &&
               string_counts[DNA_C] == seq_counts[DNA_C] &&
               string_counts[DNA_G] == seq_counts[DNA_G] &&
               string_counts[DNA_T] == seq_counts[DNA_T] &&
               string_counts[DNA_N] == seq_counts[DNA_N]
    end

    function check_nucleotide_count(::Type{RNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(RNA, seq)
        seq_counts = composition(RNASequence(seq))
        return string_counts[RNA_A] == seq_counts[RNA_A] &&
               string_counts[RNA_C] == seq_counts[RNA_C] &&
               string_counts[RNA_G] == seq_counts[RNA_G] &&
               string_counts[RNA_U] == seq_counts[RNA_U] &&
               string_counts[RNA_N] == seq_counts[RNA_N]
    end

    function check_kmer_nucleotide_count(::Type{DNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(DNA, seq)
        kmer_counts = composition(DNAKmer(seq))
        return string_counts[DNA_A] == kmer_counts[DNA_A] &&
               string_counts[DNA_C] == kmer_counts[DNA_C] &&
               string_counts[DNA_G] == kmer_counts[DNA_G] &&
               string_counts[DNA_T] == kmer_counts[DNA_T] &&
               string_counts[DNA_N] == kmer_counts[DNA_N]
    end

    function check_kmer_nucleotide_count(::Type{RNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(RNA, seq)
        kmer_counts = composition(RNAKmer(seq))
        return string_counts[RNA_A] == kmer_counts[RNA_A] &&
               string_counts[RNA_C] == kmer_counts[RNA_C] &&
               string_counts[RNA_G] == kmer_counts[RNA_G] &&
               string_counts[RNA_U] == kmer_counts[RNA_U] &&
               string_counts[RNA_N] == kmer_counts[RNA_N]
    end

    reps = 10
    for len in [1, 10, 32, 1000, 10000, 100000]
        @test all(Bool[check_nucleotide_count(DNA, random_dna(len)) for _ in 1:reps])
        @test all(Bool[check_nucleotide_count(RNA, random_rna(len)) for _ in 1:reps])
    end

    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_E] === 0
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_A] === 1
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_P] === 2
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_V] === 3

    for len in [1, 10, 32]
        @test all(Bool[check_kmer_nucleotide_count(DNA, random_dna_kmer(len)) for _ in 1:reps])
        @test all(Bool[check_kmer_nucleotide_count(RNA, random_rna_kmer(len)) for _ in 1:reps])
    end

    comp = composition(dna"ACGT")
    @test merge!(comp, composition(dna"ACGT")) === comp
    @test comp == composition(dna"ACGT"^2)

    for _ in 1:30
        m, n = rand(0:1000), rand(0:1000)
        seq1 = DNASequence(random_dna(m))
        seq2 = DNASequence(random_dna(n))
        @test composition(seq1 * seq2) == merge(composition(seq1), composition(seq2))
    end

    comp = composition(each(DNAKmer{2}, dna"ACGTACGT"))
    @test comp[DNAKmer("AC")] == 2
    @test comp[DNAKmer("CG")] == 2
    @test comp[DNAKmer("GT")] == 2
    @test comp[DNAKmer("TA")] == 1
    @test comp[DNAKmer("AA")] == 0

    comp = composition(each(DNAKmer{3}, dna"ACGTACGT"))
    @test comp[DNAKmer("ACG")] == 2
    @test comp[DNAKmer("CGT")] == 2
    @test comp[DNAKmer("GTA")] == 1
    @test comp[DNAKmer("TAC")] == 1
    @test comp[DNAKmer("AAA")] == 0
end

@testset "Kmer" begin
    reps = 10
    @testset "Construction and Conversions" begin
        @test DNACodon(DNA_A, DNA_G, DNA_T) === DNAKmer("AGT")
        @test RNACodon(RNA_A, RNA_G, RNA_U) === RNAKmer("AGU")

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

        # Check that kmers can be constructed from a BioSequence
        #   BioSequence → Kmer → BioSequence
        function check_biosequence_construction(seq::BioSequence)
            return convert(BioSequence, convert(Kmer, seq)) == seq
        end

        # Check that kmers can be constructed from an array of nucleotides
        #   Vector{NucleicAcid} → Kmer → Vector{NucleicAcid}
        function check_nucarray_kmer{T <: NucleicAcid}(seq::Vector{T})
            return convert(AbstractString, [convert(Char, c) for c in seq]) ==
                   convert(AbstractString, Kmer(seq...))
        end

        # Check that kmers in strings survive round trip conversion:
        #   String → BioSequence → Kmer → BioSequence → String
        function check_roundabout_construction(A, seq::AbstractString)
            T = eltype(A)
            return convert(AbstractString,
                       convert(BioSequence{A},
                           convert(Kmer,
                               convert(BioSequence{A}, seq)))) == uppercase(seq)
        end

        for len in [0, 1, 16, 32]
            # UInt64 conversions
            @test all(Bool[check_uint64_convertion(DNA, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])
            @test all(Bool[check_uint64_convertion(RNA, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])

            # String construction
            @test all(Bool[check_string_construction(DNA, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(RNA, random_rna_kmer(len)) for _ in 1:reps])

            # DNA/RNASequence Constructions
            @test all(Bool[check_dnasequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # BioSequence Construction
            @test all(Bool[check_biosequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # Construction from nucleotide arrays
            if len > 0
                @test all(Bool[check_nucarray_kmer(random_dna_kmer_nucleotides(len)) for _ in 1:reps])
                @test all(Bool[check_nucarray_kmer(random_rna_kmer_nucleotides(len)) for _ in 1:reps])
            end

            # Roundabout conversions
            @test all(Bool[check_roundabout_construction(DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
        end

        @test_throws Exception Kmer() # can't construct 0-mer using `Kmer()`
        @test_throws Exception Kmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
        @test_throws Exception Kmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
        @test_throws Exception Kmer(rna"ACGNU")# no Ns in kmers
        @test_throws Exception RNAKmer(rna"ACGNU")# no Ns in kmers
        @test_throws Exception Kmer(dna"ACGNT") # no Ns in kmers
        @test_throws Exception DNAKmer(dna"ACGNT") # no Ns in kmers
        @test_throws Exception Kmer(RNA_A, DNA_A) # no mixing of RNA and DNA
        @test_throws Exception Kmer(random_rna(33)) # no kmer larger than 32nt
        @test_throws Exception Kmer(random_dna(33)) # no kmer larger than 32nt
        @test_throws Exception Kmer(
                          RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 32nt
                          RNA_A, RNA_C, RNA_G, RNA_U,
                          RNA_A, RNA_C, RNA_G, RNA_U,
                          RNA_A, RNA_C, RNA_G, RNA_U,
                          RNA_A, RNA_C, RNA_G, RNA_U,
                          RNA_A, RNA_C, RNA_G, RNA_U,
                          RNA_A, RNA_C, RNA_G, RNA_U,
                          RNA_A, RNA_C, RNA_G, RNA_U,
                          RNA_A, RNA_C, RNA_G, RNA_U)
        @test_throws Exception Kmer(
                          DNA_A, DNA_C, DNA_G, DNA_T, # no kmer larger than 32nt
                          DNA_A, DNA_C, DNA_G, DNA_T,
                          DNA_A, DNA_C, DNA_G, DNA_T,
                          DNA_A, DNA_C, DNA_G, DNA_T,
                          DNA_A, DNA_C, DNA_G, DNA_T,
                          DNA_A, DNA_C, DNA_G, DNA_T,
                          DNA_A, DNA_C, DNA_G, DNA_T,
                          DNA_A, DNA_C, DNA_G, DNA_T,
                          DNA_A, DNA_C, DNA_G, DNA_T)

        @testset "From strings" begin
            @test DNAKmer("ACTG") == convert(Kmer, DNASequence("ACTG"))
            @test RNAKmer("ACUG") == convert(Kmer, RNASequence("ACUG"))

            # N is not allowed in Kmers
            @test_throws Exception DNAKmer("ACGTNACGT")
            @test_throws Exception RNAKmer("ACGUNACGU")

            # Test string literals
            @test kmer"ACTG" == convert(Kmer, DNASequence("ACTG"))
            @test isa(kmer"ACGT", DNAKmer{4})
        end
    end

    @testset "Comparisons" begin
        @testset "Equality" begin
            function check_seq_kmer_equality(len)
                a = DNAKmer(random_dna_kmer(len))
                b = convert(DNASequence, a)
                return a == b && b == a
            end

            for len in [1, 10, 32]
                @test all(Bool[check_seq_kmer_equality(len) for _ in 1:reps])
            end

            # True negatives
            @test DNAKmer("ACG") != RNAKmer("ACG")
            @test DNAKmer("T")   != RNAKmer("U")
            @test DNAKmer("AC")  != DNAKmer("AG")
            @test RNAKmer("AC")  != RNAKmer("AG")

            @test DNAKmer("ACG") != rna"ACG"
            @test DNAKmer("T")   != rna"U"
            @test DNAKmer("AC")  != dna"AG"
            @test RNAKmer("AC")  != rna"AG"

            @test rna"ACG" != DNAKmer("ACG")
            @test rna"U"   != DNAKmer("T")
            @test dna"AG"  != DNAKmer("AC")
            @test rna"AG"  != RNAKmer("AC")
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
            kmers = map(DNAKmer, ["AAAA", "AACT", "ACGT", "TGCA"])
            for x in kmers, y in kmers
                @test (x == y) == (hash(x) == hash(y))
            end
            kmers = map(RNAKmer, ["AAAA", "AACU", "ACGU", "UGCA"])
            for x in kmers, y in kmers
                @test (x == y) == (hash(x) == hash(y))
            end
        end
    end

    @testset "Length" begin
        for len in [0, 1, 16, 32]
            @test length(DNAKmer(random_dna_kmer(len))) == len
            @test length(RNAKmer(random_rna_kmer(len))) == len
        end
    end

    @testset "Arithmetic" begin
        x = DNAKmer("AA")
        @test x - 1 == x + (-1) == x - 0x01 == DNAKmer("TT")
        @test x + 1 == x - (-1) == x + 0x01 == DNAKmer("AC")

        x = DNAKmer("TT")
        @test x - 1 == x + (-1) == x - 0x01 == DNAKmer("TG")
        @test x + 1 == x - (-1) == x + 0x01 == DNAKmer("AA")

        base = DNAKmer("AAA")
        offset = 0
        nucs = "ACGT"
        for a in nucs, b in nucs, c in nucs
            @test base + offset == DNAKmer(string(a, b, c))
            offset += 1
        end
    end

    @testset "Order" begin
        @test DNAKmer("AA") < DNAKmer("AC") < DNAKmer("AG") < DNAKmer("AT") < DNAKmer("CA")
        @test RNAKmer("AA") < RNAKmer("AC") < RNAKmer("AG") < RNAKmer("AU") < RNAKmer("CA")
    end

    @testset "Access and Iterations" begin
        dna_kmer = DNAKmer("ACTG")
        rna_kmer = RNAKmer("ACUG")

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
            @test start(DNAKmer("ACTG")) == 1
            @test start(DNAKmer(""))     == 1

            @test next(DNAKmer("ACTG"), 1) == (DNA_A, 2)
            @test next(DNAKmer("ACTG"), 4) == (DNA_G, 5)

            @test  done(DNAKmer(""), 1)
            @test !done(DNAKmer("ACTG"), 1)
            @test !done(DNAKmer("ACTG"), 4)
            @test  done(DNAKmer("ACTG"), 5)
            @test !done(DNAKmer("ACTG"), -1)


            dna_vec = [DNA_A, DNA_C, DNA_T, DNA_G]
            @test all([nt === dna_vec[i] for (i, nt) in enumerate(dna_kmer)])
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
            @test start(RNAKmer("ACUG")) == 1
            @test start(RNAKmer(""))     == 1

            @test next(RNAKmer("ACUG"), 1) == (RNA_A, 2)
            @test next(RNAKmer("ACUG"), 4) == (RNA_G, 5)

            @test  done(RNAKmer(""), 1)
            @test !done(RNAKmer("ACUG"), 1)
            @test !done(RNAKmer("ACUG"), 4)
            @test  done(RNAKmer("ACUG"), 5)
            @test !done(RNAKmer("ACUG"), -1)

            rna_vec = [RNA_A, RNA_C, RNA_U, RNA_G]
            @test all([nt === rna_vec[i] for (i, nt) in enumerate(rna_kmer)])
        end
    end

    @testset "Random" begin
        @testset for k in 0:32
            for _ in 1:10
                kmer = rand(DNAKmer{k})
                @test isa(kmer, DNAKmer{k})
                kmer = rand(RNAKmer{k})
                @test isa(kmer, RNAKmer{k})
            end

            for size in [0, 1, 2, 5, 10, 100]
                @test length(rand(DNAKmer{k}, size)) == size
                @test length(rand(RNAKmer{k}, size)) == size
            end

            kmers = rand(DNAKmer{k}, 10_000)
            for i in 1:k
                a = sum([kmer[i] for kmer in kmers] .== DNA_A)
                c = sum([kmer[i] for kmer in kmers] .== DNA_C)
                g = sum([kmer[i] for kmer in kmers] .== DNA_G)
                t = sum([kmer[i] for kmer in kmers] .== DNA_T)
                @test 2200 ≤ a ≤ 2800
                @test 2200 ≤ c ≤ 2800
                @test 2200 ≤ g ≤ 2800
                @test 2200 ≤ t ≤ 2800
                @test a + c + g + t == 10_000
            end
        end
    end

    @testset "Find" begin
        kmer = DNAKmer("ACGAG")

        @test findnext(kmer, DNA_A, 1) == 1
        @test findnext(kmer, DNA_C, 1) == 2
        @test findnext(kmer, DNA_G, 1) == 3
        @test findnext(kmer, DNA_T, 1) == 0
        @test findnext(kmer, DNA_A, 2) == 4

        @test_throws BoundsError findnext(kmer, DNA_A, 0)
        @test_throws BoundsError findnext(kmer, DNA_A, 6)

        @test findprev(kmer, DNA_A, 5) == 4
        @test findprev(kmer, DNA_C, 5) == 2
        @test findprev(kmer, DNA_G, 5) == 5
        @test findprev(kmer, DNA_T, 5) == 0
        @test findprev(kmer, DNA_G, 4) == 3

        @test_throws BoundsError findprev(kmer, DNA_A, 0)
        @test_throws BoundsError findprev(kmer, DNA_A, 6)

        @test findfirst(kmer, DNA_A) == 1
        @test findfirst(kmer, DNA_G) == 3
        @test findlast(kmer, DNA_A) == 4
        @test findlast(kmer, DNA_G) == 5
    end

    @testset "Print" begin
        buf = IOBuffer()
        print(buf, DNAKmer(""))
        @test takebuf_string(buf) == ""

        buf = IOBuffer()
        print(buf, DNAKmer("ACGT"))
        @test takebuf_string(buf) == "ACGT"

        buf = IOBuffer()
        print(buf, RNAKmer("ACGU"))
        @test takebuf_string(buf) == "ACGU"
    end

    @testset "Transformations" begin
        function test_reverse(T, seq)
            revseq = reverse(Kmer{T,length(seq)}(seq))
            @test convert(AbstractString, revseq) == reverse(seq)
        end

        function test_dna_complement(seq)
            comp = complement(DNAKmer{length(seq)}(seq))
            @test convert(AbstractString, comp) == dna_complement(seq)
        end

        function test_rna_complement(seq)
            comp = complement(RNAKmer{length(seq)}(seq))
            @test convert(AbstractString, comp) == rna_complement(seq)
        end

        function test_dna_revcomp(seq)
            revcomp = reverse_complement(DNAKmer{length(seq)}(seq))
            @test convert(AbstractString, revcomp) == reverse(dna_complement(seq))
        end

        function test_rna_revcomp(seq)
            revcomp = reverse_complement(RNAKmer{length(seq)}(seq))
            @test convert(AbstractString, revcomp) == reverse(rna_complement(seq))
        end

        @testset "Reverse" begin
            for len in 0:32, _ in 1:10
                test_reverse(DNA, random_dna_kmer(len))
                test_reverse(RNA, random_rna_kmer(len))
            end

            seq = dna"AAAAAAAAAAAAAAAAAAAAAAAAAAAAGATAC"
            @test reverse(seq[(length(seq)-9):length(seq)]) == dna"CATAGAAAAA"
        end

        @testset "Complement" begin
            for len in 0:32, _ in 1:10
                test_dna_complement(random_dna_kmer(len))
                test_rna_complement(random_rna_kmer(len))
            end
        end

        @testset "Reverse Complement" begin
            for len in 0:32, _ in 1:10
                test_dna_revcomp(random_dna_kmer(len))
                test_rna_revcomp(random_rna_kmer(len))
            end
        end
    end

    @testset "Mismatches" begin
        function test_mismatches(a, b)
            count = 0
            for (x, y) in zip(a, b)
                count += x != y
            end
            @test mismatches(a, b) === mismatches(b, a) === count
        end

        for len in 0:32, _ in 1:10
            a = random_dna_kmer(len)
            b = random_dna_kmer(len)
            test_mismatches(DNAKmer(a), DNAKmer(b))

            a = random_rna_kmer(len)
            b = random_rna_kmer(len)
            test_mismatches(RNAKmer(a), RNAKmer(b))
        end
    end

    @testset "Each k-mer" begin
        function string_eachkmer(seq::AbstractString, k, step)
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

        function test_eachkmer(S, seq::AbstractString, k, step)
            xs = [convert(AbstractString, x)
                  for (i, x) in collect(each(Kmer{eltype(S),k}, S(seq), step))]
            ys = [convert(AbstractString, x)
                  for (i, x) in collect(eachkmer(S(seq), k, step))]
            zs = string_eachkmer(seq, k, step)
            @test xs == ys == zs
        end

        for k in [0, 1, 3, 16, 32], step in 1:3, len in [1, 2, 3, 5, 10, 100, 1000]
            test_eachkmer(BioSequence{DNAAlphabet{4}}, random_dna(len), k, step)
            test_eachkmer(BioSequence{RNAAlphabet{4}}, random_rna(len), k, step)
            test_eachkmer(ReferenceSequence, random_dna(len), k, step)

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_eachkmer(BioSequence{DNAAlphabet{2}}, random_dna(len, probs), k, step)
            test_eachkmer(BioSequence{RNAAlphabet{2}}, random_rna(len, probs), k, step)
        end

        @test isempty(collect(each(DNAKmer{1}, dna"")))
        @test isempty(collect(each(DNAKmer{1}, dna"NNNNNNNNNN")))
        @test_throws Exception each(DNAKmer{-1}, dna"ACGT")
        @test_throws Exception each(DNAKmer{33}, dna"ACGT")
    end

    @testset "De Bruijn Neighbors" begin
        @test collect(neighbors(DNAKmer("ACG")))  ==
            DNAKmer["CGA",  "CGC",  "CGG",  "CGT" ]
        @test collect(neighbors(DNAKmer("GGGG"))) ==
            DNAKmer["GGGA", "GGGC", "GGGG", "GGGT"]
        @test collect(neighbors(RNAKmer("ACG")))  ==
            RNAKmer["CGA",  "CGC",  "CGG",  "CGU" ]
        @test collect(neighbors(RNAKmer("GGGG"))) ==
            RNAKmer["GGGA", "GGGC", "GGGG", "GGGU"]
    end

    @testset "Shuffle" begin
        for s in ["", "A", "C", "G", "T"]
            kmer = DNAKmer(s)
            @test kmer === shuffle(kmer)
        end

        function count(kmer)
            a = c = g = t = 0
            for x in kmer
                a += x == DNA_A
                c += x == DNA_C
                g += x == DNA_G
                t += x == DNA_T
            end
            return a, c, g, t
        end

        for k in 0:32, _ in 1:10
            kmer = rand(DNAKmer{k})
            @test count(kmer) == count(shuffle(kmer))
            if k ≥ 30
                @test kmer != shuffle(kmer)
            end
        end
    end
end

@testset "Search" begin
    @testset "Exact" begin
        seq = dna"ACGTACG"

        @testset "forward" begin
            @test search(seq, DNA_A) === 1:1
            @test search(seq, DNA_N) === 1:1
            @test search(seq, DNA_T) === 4:4
            @test search(seq, DNA_T, 5) === 0:-1
            @test search(seq, DNA_A, 2, 6) === 5:5

            @test search(seq, dna"") === 1:0
            @test search(seq, dna"AC") === 1:2
            @test search(seq, dna"AC", 2) === 5:6
            @test search(seq, dna"AC", 2, 5) === 0:-1
            @test search(seq, dna"TG") === 0:-1
            @test search(seq, dna"TN") === 4:5
            @test search(seq, dna"ACG") === 1:3
            @test search(seq, dna"ACG", 2) === 5:7
            @test search(seq, seq) === 1:endof(seq)

            @test search(dna"", dna"") === 1:0
            @test search(dna"", dna"", -1) === 1:0
            @test search(dna"", dna"", 2) === 0:-1

            @test searchindex(seq, dna"") === 1
            @test searchindex(seq, dna"AC") === 1
            @test searchindex(seq, dna"AC", 2) === 5
            @test searchindex(seq, dna"AC", 2, 5) === 0

            query = ExactSearchQuery(dna"ACG")
            @test search(seq, query) === 1:3
            @test search(seq, query, 2) === 5:7
            @test search(seq, query, 2, 6) === 0:-1
            @test searchindex(seq, query) === 1
            @test searchindex(seq, query, 2) === 5
            @test searchindex(seq, query, 2, 6) === 0
        end

        @testset "backward" begin
            @test rsearch(seq, DNA_A) === 5:5
            @test rsearch(seq, DNA_N) === 7:7
            @test rsearch(seq, DNA_T) === 4:4
            @test rsearch(seq, DNA_T, 3) === 0:-1
            @test rsearch(seq, DNA_T, 6, 3) === 4:4

            @test rsearch(seq, dna"") === 8:7
            @test rsearch(seq, dna"AC") === 5:6
            @test rsearch(seq, dna"AC", 5) === 1:2
            @test rsearch(seq, dna"AC", 5, 2) === 0:-1
            @test rsearch(seq, dna"TG") === 0:-1
            @test rsearch(seq, dna"TN") === 4:5
            @test rsearch(seq, dna"ACG") === 5:7
            @test rsearch(seq, dna"ACG", 6) === 1:3
            @test rsearch(seq, seq) === 1:endof(seq)

            @test rsearch(dna"", dna"") === 1:0
            @test rsearch(dna"", dna"", 2) === 1:0
            @test rsearch(dna"", dna"", -1) === 0:-1

            @test rsearchindex(seq, dna"") === 8
            @test rsearchindex(seq, dna"AC") === 5
            @test rsearchindex(seq, dna"AC", 5) === 1
            @test rsearchindex(seq, dna"AC", 5, 2) === 0

            query = ExactSearchQuery(dna"ACG")
            @test rsearch(seq, query) === 5:7
            @test rsearch(seq, query, 6) === 1:3
            @test rsearch(seq, query, 2, 6) === 0:-1
            @test rsearchindex(seq, query) === 5
            @test rsearchindex(seq, query, 6) === 1
            @test rsearchindex(seq, query, 6, 2) === 0
        end
    end

    @testset "Approximate" begin
        seq = dna"ACGTACG"

        @testset "forward" begin
            @test approxsearch(seq, dna"", 0) === 1:0
            @test approxsearch(seq, dna"AC", 0) === 1:2
            @test approxsearch(seq, dna"AC", 0, 2) === 5:6
            @test approxsearch(seq, dna"AC", 0, 2, 5) === 0:-1
            @test approxsearch(seq, dna"AT", 0) === 0:-1
            @test approxsearch(seq, dna"AT", 1) === 1:1
            @test approxsearch(seq, dna"AT", 1, 2) === 3:4
            @test approxsearch(seq, dna"AT", 1, 2, 3) === 0:-1
            @test approxsearch(seq, dna"NG", 0) === 2:3
            @test approxsearch(seq, dna"NG", 1) === 1:1
            @test approxsearch(seq, dna"GN", 0) === 3:4
            @test approxsearch(seq, dna"GN", 1) === 1:1
            @test approxsearch(seq, dna"ACG", 0) === 1:3
            @test approxsearch(seq, dna"ACG", 1) === 1:2
            @test approxsearch(seq, dna"ACG", 2) === 1:1
            @test approxsearch(seq, dna"ACG", 3) === 1:0
            @test approxsearch(seq, dna"ACG", 4) === 1:0

            @test approxsearchindex(seq, dna"", 0) === 1
            @test approxsearchindex(seq, dna"AC", 0) === 1
            @test approxsearchindex(seq, dna"AC", 0, 2) === 5
            @test approxsearchindex(seq, dna"AC", 0, 2, 5) === 0

            query = ApproximateSearchQuery(dna"ACG")
            @test approxsearch(seq, query, 1) === 1:2
            @test approxsearch(seq, query, 1, 2) === 2:3
            @test approxsearch(seq, query, 1, 2, 2) === 0:-1
            @test approxsearchindex(seq, query, 1) === 1
            @test approxsearchindex(seq, query, 1, 2) === 2
            @test approxsearchindex(seq, query, 1, 2, 2) === 0
        end

        @testset "backward" begin
            # TODO: maybe this should return 8:7 like rsearch
            @test approxrsearch(seq, dna"", 0) === 7:6
            @test approxrsearch(seq, dna"AC", 0) === 5:6
            @test approxrsearch(seq, dna"AC", 0, 5) === 1:2
            @test approxrsearch(seq, dna"AC", 0, 5, 2) === 0:-1
            @test approxrsearch(seq, dna"AT", 0) === 0:-1
            @test approxrsearch(seq, dna"AT", 1) === 5:6
            @test approxrsearch(seq, dna"AT", 1, 6) === 5:6
            @test approxrsearch(seq, dna"AT", 1, 5) === 5:5
            @test approxrsearch(seq, dna"AT", 1, 3, 2) === 0:-1
            @test approxrsearch(seq, dna"NG", 0) === 6:7
            @test approxrsearch(seq, dna"NG", 0, 6) === 2:3
            @test approxrsearch(seq, dna"GN", 0) === 3:4
            @test approxrsearch(seq, dna"GN", 1) === 7:7
            @test approxrsearch(seq, dna"ACG", 0) === 5:7
            @test approxrsearch(seq, dna"ACG", 1) === 6:7
            @test approxrsearch(seq, dna"ACG", 2) === 7:7
            @test approxrsearch(seq, dna"ACG", 3) === 7:6
            @test approxrsearch(seq, dna"ACG", 4) === 7:6

            # TODO: maybe this should return 8 like rsearchindex
            @test approxrsearchindex(seq, dna"", 0) === 7
            @test approxrsearchindex(seq, dna"AC", 0) === 5
            @test approxrsearchindex(seq, dna"AC", 0, 5) === 1
            @test approxrsearchindex(seq, dna"AC", 0, 5, 2) === 0

            query = ApproximateSearchQuery(dna"ACG")
            @test approxrsearch(seq, query, 1, 7) === 6:7
            @test approxrsearch(seq, query, 1, 6) === 5:6
            @test approxrsearch(seq, query, 1, 6, 6) === 0:-1
            @test approxrsearchindex(seq, query, 1, 7) === 6
            @test approxrsearchindex(seq, query, 1, 6) === 5
            @test approxrsearchindex(seq, query, 1, 6, 6) === 0
        end
    end

    @testset "Regular Expression" begin
        @test isa(biore"A+"d, Seq.RE.Regex{DNA})
        @test isa(biore"A+"r, Seq.RE.Regex{RNA})
        @test isa(biore"A+"a, Seq.RE.Regex{AminoAcid})
        @test isa(biore"A+"dna, Seq.RE.Regex{DNA})
        @test isa(biore"A+"rna, Seq.RE.Regex{RNA})
        @test isa(biore"A+"aa, Seq.RE.Regex{AminoAcid})
        @test_throws Exception eval(:(biore"A+"))
        @test_throws Exception eval(:(biore"A+"foo))
        @test string(biore"A+"dna) == "biore\"A+\"dna"
        @test string(biore"A+"rna) == "biore\"A+\"rna"
        @test string(biore"A+"aa) == "biore\"A+\"aa"

        @test  ismatch(biore"A"d, dna"A")
        @test !ismatch(biore"A"d, dna"C")
        @test !ismatch(biore"A"d, dna"G")
        @test !ismatch(biore"A"d, dna"T")
        @test  ismatch(biore"N"d, dna"A")
        @test  ismatch(biore"N"d, dna"C")
        @test  ismatch(biore"N"d, dna"G")
        @test  ismatch(biore"N"d, dna"T")
        @test  ismatch(biore"[AT]"d, dna"A")
        @test !ismatch(biore"[AT]"d, dna"C")
        @test !ismatch(biore"[AT]"d, dna"G")
        @test  ismatch(biore"[AT]"d, dna"T")
        @test !ismatch(biore"[^AT]"d, dna"A")
        @test  ismatch(biore"[^AT]"d, dna"C")
        @test  ismatch(biore"[^AT]"d, dna"G")
        @test !ismatch(biore"[^AT]"d, dna"T")

        re = biore"^A(C+G*)(T{2,})N$"d
        @test !ismatch(re, dna"AC")
        @test !ismatch(re, dna"AGTT")
        @test !ismatch(re, dna"CCGTT")
        @test !ismatch(re, dna"ACTT")
        @test !ismatch(re, dna"ACTTGT")
        @test  ismatch(re, dna"ACGTTA")
        @test  ismatch(re, dna"ACGTTT")
        @test  ismatch(re, dna"ACCGGTTT")
        @test  ismatch(re, dna"ACCGGTTT")
        @test  ismatch(re, dna"ACCGGTTTA")
        @test  ismatch(re, dna"ACCGGTTTG")

        @test matched(match(re, dna"ACCGTTTTA")) == dna"ACCGTTTTA"
        @test get(captured(match(re, dna"ACCGTTTTA"))[1]) == dna"CCG"
        @test get(captured(match(re, dna"ACCGTTTTA"))[2]) == dna"TTTT"

        # greedy
        @test matched(match(biore"A*"d, dna"AAA")) == dna"AAA"
        @test matched(match(biore"A+"d, dna"AAA")) == dna"AAA"
        @test matched(match(biore"A?"d, dna"AAA")) == dna"A"
        @test matched(match(biore"A{2}"d, dna"AAA")) == dna"AA"
        @test matched(match(biore"A{2,}"d, dna"AAA")) == dna"AAA"
        @test matched(match(biore"A{2,4}"d, dna"AAA")) == dna"AAA"

        # lazy
        @test matched(match(biore"A*?"d, dna"AAA")) == dna""
        @test matched(match(biore"A+?"d, dna"AAA")) == dna"A"
        @test matched(match(biore"A??"d, dna"AAA")) == dna""
        @test matched(match(biore"A{2}?"d, dna"AAA")) == dna"AA"
        @test matched(match(biore"A{2,}?"d, dna"AAA")) == dna"AA"
        @test matched(match(biore"A{2,4}?"d, dna"AAA")) == dna"AA"

        # search
        @test search(dna"ACGTAAT", biore"A+"d) == 1:1
        @test search(dna"ACGTAAT", biore"A+"d, 1) == 1:1
        @test search(dna"ACGTAAT", biore"A+"d, 2) == 5:6
        @test search(dna"ACGTAAT", biore"A+"d, 7) == 0:-1

        # eachmatch
        matches = [dna"CG", dna"GC", dna"GC", dna"CG"]
        for (i, m) in enumerate(eachmatch(biore"GC|CG"d, dna"ACGTTATGCATGGCG"))
            @test matched(m) == matches[i]
        end
        matches = [dna"CG", dna"GC", dna"GC"]
        for (i, m) in enumerate(collect(eachmatch(biore"GC|CG"d, dna"ACGTTATGCATGGCG", false)))
            @test matched(m) == matches[i]
        end

        # matchall
        @test matchall(biore"A*"d, dna"") == [dna""]
        @test matchall(biore"A*"d, dna"AAA") == [
            dna"AAA", dna"AA", dna"A", dna"",
            dna"AA",  dna"A",  dna"",
            dna"A",   dna""]
        @test matchall(biore"AC*G*T"d, dna"ACCGGGT") == [dna"ACCGGGT"]

        @test matchall(biore"A*"d, dna"", false) == [dna""]
        @test matchall(biore"A*"d, dna"AAA", false) == [dna"AAA"]
        @test matchall(biore"AC*G*T"d, dna"ACCGGGT", false) == [dna"ACCGGGT"]

        # RNA and Amino acid
        @test  ismatch(biore"U(A[AG]|GA)$"r, rna"AUUGUAUGA")
        @test !ismatch(biore"U(A[AG]|GA)$"r, rna"AUUGUAUGG")
        @test  ismatch(biore"T+[NQ]A?P"a, aa"MTTQAPMFTQPL")
        @test  ismatch(biore"T+[NQ]A?P"a, aa"MTTAAPMFTQPL")
        @test !ismatch(biore"T+[NQ]A?P"a, aa"MTTAAPMFSQPL")

        # PROSITE
        @test  ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"ADVAARRK")
        @test  ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVAARRK")
        @test !ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"ADVAARRE")
        @test !ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVAARK")
        @test  ismatch(prosite"<[AC]-x-V-x(4)-{ED}>", aa"ADVAARRK")
        @test !ismatch(prosite"<[AC]-x-V-x(4)-{ED}>", aa"AADVAARRK")
        @test !ismatch(prosite"<[AC]-x-V-x(4)-{ED}>", aa"ADVAARRKA")
    end
end

@testset "Translation" begin
    # crummy string translation to test against
    standard_genetic_code_dict = Dict{String,Char}(
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
    for len in [1, 10, 32, 1000, 10000, 100000]
        @test all(Bool[check_translate(random_translatable_rna(len)) for _ in 1:reps])
    end

    # ambiguous codons
    @test translate(rna"YUGMGG") == aa"LR"
    @test translate(rna"GAYGARGAM") == aa"DEX"

    @test_throws Exception translate(dna"ACGTACGTA") # can't translate DNA
    @test_throws Exception translate(rna"ACGUACGU")  # can't translate non-multiples of three
    # can't translate N
    @test_throws Exception translate(rna"ACGUACGNU", allow_ambiguous_codons=false)

    # issue #133
    @test translate(rna"GAN") == aa"X"
end

@testset "Demultiplexer" begin
    function randdna(n)
        return DNASequence(rand([DNA_A, DNA_C, DNA_G, DNA_T, DNA_N], n))
    end

    function make_errors(seq, p=0.03)
        seq = copy(seq)
        nucs = DNA['A', 'C', 'G', 'T', 'N']
        i = 1
        while i ≤ endof(seq)
            if rand() < p
                r = rand()
                if r < 0.1
                    # insertion
                    insert!(seq, i, rand(nucs))
                elseif r < 0.2
                    # deletion
                    deleteat!(seq, i)
                else
                    # substitution
                    seq[i] = rand(nucs)
                end
            end
            i += 1
        end
        return seq
    end

    @testset "Hamming distance" begin
        barcodes = DNASequence["ATGG", "CAGA", "GGAA", "TACG"]
        dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming)

        for i in 1:endof(barcodes)
            @test dplxr[i] == barcodes[i]
        end

        @test demultiplex(dplxr, dna"ATGG") === (1, 0)
        @test demultiplex(dplxr, dna"CAGA") === (2, 0)
        @test demultiplex(dplxr, dna"GGAA") === (3, 0)
        @test demultiplex(dplxr, dna"TACG") === (4, 0)

        # every 1bp substitution is recoverable
        for (i, barcode) in enumerate(barcodes)
            # substitution
            for j in 1:endof(barcode)
                barcode′ = copy(barcode)
                for nt in [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
                    barcode′[j] = nt
                    @test demultiplex(dplxr, barcode′) == (i, mismatches(barcode, barcode′))
                end
            end
        end

        n_ok = 0
        n_ok_with_fallback = 0
        for _ in 1:10_000
            i = rand(1:endof(barcodes))
            seq = make_errors(barcodes[i] * randdna(10))
            n_ok += demultiplex(dplxr, seq)[1] == i
            n_ok_with_fallback += demultiplex(dplxr, seq, true)[1] == i
        end
        # empirically, n_ok / 10_000 is ~0.985
        # @show n_ok
        @test n_ok / 10_000 > 0.98
        @test n_ok < n_ok_with_fallback
    end

    @testset "Levenshtein distance" begin
        barcodes = DNASequence["ATGG", "CAGA", "GGAA", "TACG"]
        dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:levenshtein)

        for i in 1:endof(barcodes)
            @test dplxr[i] == barcodes[i]
        end

        @test demultiplex(dplxr, dna"ATGG") === (1, 0)
        @test demultiplex(dplxr, dna"CAGA") === (2, 0)
        @test demultiplex(dplxr, dna"GGAA") === (3, 0)
        @test demultiplex(dplxr, dna"TACG") === (4, 0)

        # every 1bp substitution/insertion/deletion is recoverable
        for (i, barcode) in enumerate(barcodes)
            # substitution
            for j in 1:endof(barcode)
                barcode′ = copy(barcode)
                for nt in [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
                    barcode′[j] = nt
                    i′, d = demultiplex(dplxr, barcode′)
                    @test i′ == i && 0 ≤ d ≤ 1
                end
            end
            # insertion
            for j in 1:endof(barcode)
                barcode′ = copy(barcode)
                insert!(barcode′, j, rand([DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]))
                i′, d = demultiplex(dplxr, barcode′)
                @test i′ == i && 0 ≤ d ≤ 1
            end
            # deletion
            for j in 1:endof(barcode)
                barcode′ = copy(barcode)
                deleteat!(barcode′, j)
                i′, d = demultiplex(dplxr, barcode′)
                @test i′ == i && 0 ≤ d ≤ 1
            end
        end

        n_ok = 0
        n_ok_with_fallback = 0
        for _ in 1:10_000
            i = rand(1:endof(barcodes))
            seq = make_errors(barcodes[i] * randdna(10))
            n_ok += demultiplex(dplxr, seq)[1] == i
            n_ok_with_fallback += demultiplex(dplxr, seq, true)[1] == i
        end
        # empirically, n_ok / 10_000 is ~0.995
        # @show n_ok
        @test n_ok / 10_000 > 0.99
        @test n_ok < n_ok_with_fallback
    end
end

@testset "SeqRecord" begin
    rec = SeqRecord("seq1", dna"ACGTN")
    @test isa(rec, SeqRecord{DNASequence,Void})
    @test seqname(rec) == "seq1"
    @test sequence(rec) == dna"ACGTN"
    @test metadata(rec) == nothing

    buf = IOBuffer()
    show(buf, rec)
    @test startswith(takebuf_string(buf), "Bio.Seq.SeqRecord")

    rec1 = SeqRecord("seq1", dna"ACGTN")
    rec2 = SeqRecord("seq2", dna"ACGTN")
    rec3 = SeqRecord("seq1", dna"GGAGT")
    @test rec == rec1
    @test rec != rec2
    @test rec != rec3
    @test rec == copy(rec)
end

@testset "Reading and Writing" begin
    @testset "FASTA" begin
        @testset "Record" begin
            record = FASTA.Record()
            @test !isfilled(record)
            @test_throws ArgumentError FASTA.identifier(record)

            record = FASTA.Record(b">foo\nACGT\n")
            @test isfilled(record)
            @test hasseqname(record)
            @test FASTA.hasidentifier(record)
            @test seqname(record) == FASTA.identifier(record) == "foo"
            @test !FASTA.hasdescription(record)
            @test_throws MissingFieldException FASTA.description(record)
            @test hassequence(record)
            @test FASTA.hassequence(record)
            @test FASTA.sequence(record) == dna"ACGT"
            @test FASTA.sequence(record, 2:3) == dna"CG"
            @test FASTA.sequence(String, record) == "ACGT"
            @test FASTA.sequence(String, record, 2:3) == "CG"
            @test record == FASTA.Record(">foo\nACGT\n")

            record = FASTA.Record(b"""
            >CYS1_DICDI fragment
            SCWSFSTTGNVEGQHFISQNKL
            VSLSEQNLVDCDHECMEYEGE
            """)
            @test isfilled(record)
            @test FASTA.identifier(record) == "CYS1_DICDI"
            @test FASTA.description(record) == "fragment"
            @test FASTA.sequence(record) == aa"SCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGE"
            @test FASTA.sequence(record, 10:15) == aa"NVEGQH"
        end

        output = IOBuffer()
        writer = FASTA.Writer(output, 5)
        write(writer, FASTA.Record("seq1", dna"TTA"))
        write(writer, FASTA.Record("seq2", "some description", dna"ACGTNN"))
        flush(writer)
        @test takebuf_string(output) == """
        >seq1
        TTA
        >seq2 some description
        ACGTN
        N
        """

        reader = FASTA.Reader(IOBuffer(
        """
        >seqA some description
        QIKDLLVSSSTDLDTTLKMK
        ILELPFASGDLSM
        >seqB
        VLMALGMTDLFIPSANLTG*
        """))
        record = FASTA.Record()
        @test read!(reader, record) === record
        @test FASTA.identifier(record) == "seqA"
        @test FASTA.description(record) == "some description"
        @test FASTA.sequence(record) == aa"QIKDLLVSSSTDLDTTLKMKILELPFASGDLSM"
        @test read!(reader, record) === record
        @test FASTA.identifier(record) == "seqB"
        @test !FASTA.hasdescription(record)
        @test FASTA.sequence(record) == aa"VLMALGMTDLFIPSANLTG*"

        function test_fasta_parse(filename, valid)
            # Reading from a stream
            stream = open(FASTA.Reader, filename)
            @test eltype(stream) == FASTA.Record
            if valid
                for seqrec in stream end
                @test true  # no error
                @test close(stream) === nothing
            else
                @test_throws Exception begin
                    for seqrec in stream end
                end
                return
            end

            # in-place parsing
            stream = open(FASTA.Reader, filename)
            entry = eltype(stream)()
            while !eof(stream)
                read!(stream, entry)
            end

            # Check round trip
            output = IOBuffer()
            writer = FASTA.Writer(output, width=60)
            expected_entries = Any[]
            for seqrec in open(FASTA.Reader, filename)
                write(writer, seqrec)
                push!(expected_entries, seqrec)
            end
            flush(writer)

            seekstart(output)
            read_entries = FASTA.Record[]
            for seqrec in FASTA.Reader(output)
                push!(read_entries, seqrec)
            end

            @test expected_entries == read_entries
        end

        get_bio_fmt_specimens()
        path = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "FASTA")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            tags = specimen["tags"]
            valid = get(specimen, "valid", true)
            if contains(tags, "comments")
                # currently comments are not supported
                continue
            end
            test_fasta_parse(joinpath(path, specimen["filename"]), valid)
        end

        @testset "Faidx" begin
            fastastr = """
            >chr1
            CCACACCACACCCACACACC
            >chr2
            ATGCATGCATGCAT
            GCATGCATGCATGC
            >chr3
            AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGA
            TGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATA
            >chr4
            TACTT
            """
            # generated with `samtools faidx`
            faistr = """
            chr1	20	6	20	21
            chr2	28	33	14	15
            chr3	100	69	50	51
            chr4	5	177	5	6
            """
            mktempdir() do dir
                filepath = joinpath(dir, "test.fa")
                write(filepath, fastastr)
                write(filepath * ".fai", faistr)
                open(FASTA.Reader, filepath, index=filepath * ".fai") do reader
                    chr3 = reader["chr3"]
                    @test FASTA.identifier(chr3) == "chr3"
                    @test FASTA.sequence(chr3) == dna"""
                    AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGA
                    TGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATA
                    """

                    chr2 = reader["chr2"]
                    @test FASTA.identifier(chr2) == "chr2"
                    @test FASTA.sequence(chr2) == dna"""
                    ATGCATGCATGCAT
                    GCATGCATGCATGC
                    """

                    chr4 = reader["chr4"]
                    @test FASTA.identifier(chr4) == "chr4"
                    @test FASTA.sequence(chr4) == dna"""
                    TACTT
                    """

                    chr1 = reader["chr1"]
                    @test FASTA.identifier(chr1) == "chr1"
                    @test FASTA.sequence(chr1) == dna"""
                    CCACACCACACCCACACACC
                    """

                    @test_throws ArgumentError reader["chr5"]
                end
            end

            # invalid index
            @test_throws ArgumentError FASTA.Reader(IOBuffer(fastastr), index=π)
        end

        @testset "append" begin
            intempdir() do
                filepath = "test.fa"
                writer = open(FASTA.Writer, filepath)
                write(writer, FASTA.Record("seq1", dna"AAA"))
                close(writer)
                writer = open(FASTA.Writer, filepath, append=true)
                write(writer, FASTA.Record("seq2", dna"CCC"))
                close(writer)
                seqs = open(collect, FASTA.Reader, filepath)
                @test length(seqs) == 2
            end
        end
    end

    @testset "FASTQ" begin
        @test isa(FASTQ.Record("1", dna"AA", UInt8[10, 11]), FASTQ.Record)
        @test isa(FASTQ.Record("1", "desc.", dna"AA", UInt8[10, 11]), FASTQ.Record)
        @test_throws ArgumentError FASTQ.Record("1", dna"AA", UInt8[10])

        output = IOBuffer()
        writer = FASTQ.Writer(output)
        write(writer, FASTQ.Record("1", dna"AN", UInt8[11, 25]))
        write(writer, FASTQ.Record("2", "high quality", dna"TGA", UInt8[40, 41, 45]))
        flush(writer)
        @test takebuf_string(output) == """
        @1
        AN
        +
        ,:
        @2 high quality
        TGA
        +
        IJN
        """

        output = IOBuffer()
        writer = FASTQ.Writer(output, quality_header=true)
        write(writer, FASTQ.Record("1", dna"AN", UInt8[11, 25]))
        write(writer, FASTQ.Record("2", "high quality", dna"TGA", UInt8[40, 41, 45]))
        flush(writer)
        @test takebuf_string(output) == """
        @1
        AN
        +1
        ,:
        @2 high quality
        TGA
        +2 high quality
        IJN
        """

        @testset "Record" begin
            record = FASTQ.Record()
            @test !isfilled(record)

            record = FASTQ.Record("""
            @SRR1238088.1.1 HWI-ST499:111:D0G94ACXX:1:1101:1173:2105
            AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA
            +SRR1238088.1.1 HWI-ST499:111:D0G94ACXX:1:1101:1173:2105
            @BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ
            """)
            @test isfilled(record)
            @test FASTQ.hasidentifier(record) == hasseqname(record) == true
            @test FASTQ.identifier(record) == seqname(record) == "SRR1238088.1.1"
            @test FASTQ.hasdescription(record)
            @test FASTQ.description(record) == "HWI-ST499:111:D0G94ACXX:1:1101:1173:2105"
            @test FASTQ.hassequence(record) == hassequence(record) == true
            @test FASTQ.sequence(DNASequence, record) == dna"AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
            @test FASTQ.sequence(record) == sequence(record) == dna"AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
            @test FASTQ.sequence(String, record) == "AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
            @test FASTQ.hasquality(record)
            @test FASTQ.quality(record) == b"@BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ" .- 33

            record = FASTQ.Record("""
            @SRR1238088.1.1
            AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA
            +
            @BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ
            """)
            @test isfilled(record)
            @test !FASTQ.hasdescription(record)
        end

        function test_records(rs1, rs2)
            if length(rs1) != length(rs2)
                return false
            end
            for (r1, r2) in zip(rs1, rs2)
                if FASTQ.identifier(r1) != FASTQ.identifier(r2) ||
                   FASTQ.sequence(r1)   != FASTQ.sequence(r2)   ||
                   FASTQ.quality(r1)    != FASTQ.quality(r2)
                    return false
                end
            end
            return true
        end

        function test_fastq_parse(filename, valid)
            # Reading from a reader
            reader = open(FASTQ.Reader, filename)
            @test eltype(reader) == FASTQ.Record
            if valid
                for record in reader end
                @test true  # no error
                @test close(reader) === nothing
            else
                @test_throws Exception begin
                    for record in reader end
                end
                return
            end

            # in-place parsing
            reader = open(FASTQ.Reader, filename)
            record = eltype(reader)()
            try
                while true
                    read!(reader, record)
                end
            catch ex
                close(reader)
                if !isa(ex, EOFError)
                    rethrow()
                end
            end

            # Check round trip
            output = IOBuffer()
            writer = FASTQ.Writer(output)
            expected_entries = FASTQ.Record[]
            for record in open(FASTQ.Reader, filename)
                write(writer, record)
                push!(expected_entries, record)
            end
            flush(writer)

            seekstart(output)
            read_entries = FASTQ.Record[]
            for record in FASTQ.Reader(output)
                push!(read_entries, record)
            end

            return test_records(expected_entries, read_entries)
        end

        get_bio_fmt_specimens()
        path = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "FASTQ")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            tags = split(get(specimen, "tags", ""))
            valid = get(specimen, "valid", true)
            # currently unsupported features
            if any(t ∈ tags for t in ["gaps", "rna", "comments", "linewrap"])
                continue
            end
            filename = specimen["filename"]
            test_fastq_parse(joinpath(path, filename), valid)
        end

        @testset "invalid quality encoding" begin
            # Sanger full range (note escape characters before '$' and '\')
            record = FASTQ.Record("""
            @FAKE0001 Original version has PHRED scores from 0 to 93 inclusive (in that order)
            ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
            +
            !"#\$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~
            """)

            # the range is not enough in these encodings
            for encoding in (:solexa, :illumina13, :illumina15)
                @test_throws ErrorException FASTQ.quality(record, encoding)
            end

            # the range is enough in these encodings
            for encoding in (:sanger, :illumina18)
                @test FASTQ.quality(record, encoding) == collect(0:93)
            end
        end

        @testset "fill ambiguous nucleotides" begin
            input = IOBuffer("""
            @seq1
            ACGTNRacgtnr
            +
            BBBB##AAAA##
            """)
            @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=nothing))) == dna"ACGTNRACGTNR"
            seekstart(input)
            @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_A)))   == dna"ACGTAAACGTAA"
            seekstart(input)
            @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_G)))   == dna"ACGTGGACGTGG"
            seekstart(input)
            @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_N)))   == dna"ACGTNNACGTNN"
            seekstart(input)
            @test FASTQ.sequence(BioSequence{DNAAlphabet{2}}, first(FASTQ.Reader(input, fill_ambiguous=DNA_A))) == dna"ACGTAAACGTAA"
        end
    end

    @testset "2bit" begin
        buffer = IOBuffer()
        writer = TwoBit.Writer(buffer, ["chr1", "chr2"])
        chr1 = dna"ACGTNN"
        chr2 = dna"N"^100 * dna"ACGT"^100 * dna"N"^100
        write(writer, SeqRecord("chr1", chr1))
        write(writer, SeqRecord("chr2", chr2))
        seekstart(buffer)
        reader = TwoBit.Reader(buffer)
        @test length(reader) == 2
        @test TwoBit.seqnames(reader) == ["chr1", "chr2"]
        @test TwoBit.sequence(reader["chr1"]) == chr1
        @test TwoBit.sequence(reader["chr2"]) == chr2
        @test TwoBit.sequence(reader["chr1"]) == chr1
        @test_throws KeyError reader["chr10"]
        @test_throws BoundsError reader[3]

        function check_2bit_parse(filename)
            stream = open(TwoBit.Reader, filename)
            @test eltype(stream) === TwoBit.Record
            # read from a stream
            for record in stream
                @test hassequence(record) == TwoBit.hassequence(record) == true
                @test TwoBit.sequence(ReferenceSequence, record) == TwoBit.sequence(DNASequence, record)
            end
            close(stream)

            # round trip
            buffer = IOBuffer()
            reader = open(TwoBit.Reader, filename)
            writer = TwoBit.Writer(buffer, TwoBit.seqnames(reader))
            expected_entries = TwoBit.Record[]
            for (name, record) in zip(TwoBit.seqnames(reader), reader)
                write(writer, SeqRecord(name, TwoBit.sequence(record), TwoBit.maskedblocks(record)))
                push!(expected_entries, record)
            end

            read_entries = TwoBit.Record[]
            seekstart(buffer)
            for record in TwoBit.Reader(buffer)
                push!(read_entries, record)
            end

            return expected_entries == read_entries
        end

        get_bio_fmt_specimens()
        path = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "2bit")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            valid = get(specimen, "valid", true)
            filepath = joinpath(path, specimen["filename"])
            if valid
                @test check_2bit_parse(filepath)
            else
                @test_throws Exception check_fasta_parse(filepath)
            end
        end
    end

    @testset "ABIF Reader" begin
        function check_abif_parse(filename)
            stream = open(AbifReader, filename)

            for record in stream end
            for (a,b) in collect(stream[1]) end
            for (a,b) in getindex(stream, get_tags(stream)) end
            @test typeof(stream) == AbifReader{IOStream}
            @test get_tags(stream)[1].name == "AEPt"
            @test get_tags(stream, "DATA")[1].name == "DATA"
            @test length(stream["DATA"]) == 12
            @test length(stream[1]) == 1

            @test typeof(getindex(stream, get_tags(stream))) == Dict{String,Any}
            @test tagelements(stream, "DATA") == 12
        end

        get_bio_fmt_specimens()
        path = Pkg.dir("Bio", "test", "BioFmtSpecimens", "ABI")
        for specimen in YAML.load_file(joinpath(path, "index.yml"))
            valid = get(specimen, "valid", true)
            filepath = joinpath(path, specimen["filename"])
            if valid
                check_abif_parse(filepath)
            else
                @test_throws Exception check_abif_parse(filepath)
            end
        end
    end
end

@testset "Quality scores" begin
    @testset "Decoding base quality scores" begin
        function test_decode(encoding, values, expected)
            result = Array{Int8}(length(expected))
            Seq.decode_quality_string!(encoding, values, result, 1, length(result))
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

    @testset "Encoding base quality scores" begin
        function test_encode(encoding, values, expected)
            result = Array{UInt8}(length(expected))
            Seq.encode_quality_string!(encoding, values, result, 1, length(result))
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
end

@testset "majorityvote" begin
    dnas = [dna"CTCGATCGATCC", dna"CTCGAAAAATCA", dna"ATCGAAAAATCG", dna"ATCGGGGGATCG"]
    dna2 = [dna"CTCGATCGATCC", dna"CTCGAAA", dna"ATCGAAAAATCG", dna"ATCGGGGGATCG"]
    rnas = [rna"CUCGAUCGAUCC", rna"CUCGAAAAAUCA", rna"AUCGAAAAAUCG", rna"AUCGGGGGAUCG"]
    rna2 = [rna"CUCGAUCGAUCC", rna"CUCGAAAAAUCA", rna"AUCGAAAAAUCG", rna"AUCGUCG"]
    @test majorityvote(dnas) == dna"MTCGAAARATCG"
    @test majorityvote(rnas) == rna"MUCGAAARAUCG"
    @test_throws ArgumentError majorityvote(dna2)
    @test_throws ArgumentError majorityvote(rna2)
end

end # TestSeq
