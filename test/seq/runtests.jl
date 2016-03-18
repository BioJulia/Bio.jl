module TestSeq

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

import Bio
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
    return convert(ASCIIString, seqc)
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
    return convert(ASCIIString, seqc)
end

function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
end

@testset "Nucleotides" begin
    @testset "Conversions" begin
        @testset "UInt8" begin
            @testset "DNA conversions from UInt8" begin
                @test convert(DNANucleotide, UInt8( 0)) == DNA_A
                @test convert(DNANucleotide, UInt8( 1)) == DNA_C
                @test convert(DNANucleotide, UInt8( 2)) == DNA_G
                @test convert(DNANucleotide, UInt8( 3)) == DNA_T
                @test convert(DNANucleotide, UInt8( 4)) == DNA_M
                @test convert(DNANucleotide, UInt8( 5)) == DNA_R
                @test convert(DNANucleotide, UInt8( 6)) == DNA_W
                @test convert(DNANucleotide, UInt8( 7)) == DNA_S
                @test convert(DNANucleotide, UInt8( 8)) == DNA_Y
                @test convert(DNANucleotide, UInt8( 9)) == DNA_K
                @test convert(DNANucleotide, UInt8(10)) == DNA_V
                @test convert(DNANucleotide, UInt8(11)) == DNA_H
                @test convert(DNANucleotide, UInt8(12)) == DNA_D
                @test convert(DNANucleotide, UInt8(13)) == DNA_B
                @test convert(DNANucleotide, UInt8(14)) == DNA_N
                @test convert(DNANucleotide, UInt8(15)) == DNA_Gap
            end

            @testset "RNA conversions from UInt8" begin
                @test convert(RNANucleotide, UInt8( 0)) == RNA_A
                @test convert(RNANucleotide, UInt8( 1)) == RNA_C
                @test convert(RNANucleotide, UInt8( 2)) == RNA_G
                @test convert(RNANucleotide, UInt8( 3)) == RNA_U
                @test convert(RNANucleotide, UInt8( 4)) == RNA_M
                @test convert(RNANucleotide, UInt8( 5)) == RNA_R
                @test convert(RNANucleotide, UInt8( 6)) == RNA_W
                @test convert(RNANucleotide, UInt8( 7)) == RNA_S
                @test convert(RNANucleotide, UInt8( 8)) == RNA_Y
                @test convert(RNANucleotide, UInt8( 9)) == RNA_K
                @test convert(RNANucleotide, UInt8(10)) == RNA_V
                @test convert(RNANucleotide, UInt8(11)) == RNA_H
                @test convert(RNANucleotide, UInt8(12)) == RNA_D
                @test convert(RNANucleotide, UInt8(13)) == RNA_B
                @test convert(RNANucleotide, UInt8(14)) == RNA_N
                @test convert(RNANucleotide, UInt8(15)) == RNA_Gap
            end

            @testset "DNA conversions to UInt8" begin
                @test convert(UInt8, DNA_A)   == UInt8( 0)
                @test convert(UInt8, DNA_C)   == UInt8( 1)
                @test convert(UInt8, DNA_G)   == UInt8( 2)
                @test convert(UInt8, DNA_T)   == UInt8( 3)
                @test convert(UInt8, DNA_M)   == UInt8( 4)
                @test convert(UInt8, DNA_R)   == UInt8( 5)
                @test convert(UInt8, DNA_W)   == UInt8( 6)
                @test convert(UInt8, DNA_S)   == UInt8( 7)
                @test convert(UInt8, DNA_Y)   == UInt8( 8)
                @test convert(UInt8, DNA_K)   == UInt8( 9)
                @test convert(UInt8, DNA_V)   == UInt8(10)
                @test convert(UInt8, DNA_H)   == UInt8(11)
                @test convert(UInt8, DNA_D)   == UInt8(12)
                @test convert(UInt8, DNA_B)   == UInt8(13)
                @test convert(UInt8, DNA_N)   == UInt8(14)
                @test convert(UInt8, DNA_Gap) == UInt8(15)
            end

            @testset "RNA conversions to UInt8" begin
                @test convert(UInt8, RNA_A)   == UInt8( 0)
                @test convert(UInt8, RNA_C)   == UInt8( 1)
                @test convert(UInt8, RNA_G)   == UInt8( 2)
                @test convert(UInt8, RNA_U)   == UInt8( 3)
                @test convert(UInt8, RNA_M)   == UInt8( 4)
                @test convert(UInt8, RNA_R)   == UInt8( 5)
                @test convert(UInt8, RNA_W)   == UInt8( 6)
                @test convert(UInt8, RNA_S)   == UInt8( 7)
                @test convert(UInt8, RNA_Y)   == UInt8( 8)
                @test convert(UInt8, RNA_K)   == UInt8( 9)
                @test convert(UInt8, RNA_V)   == UInt8(10)
                @test convert(UInt8, RNA_H)   == UInt8(11)
                @test convert(UInt8, RNA_D)   == UInt8(12)
                @test convert(UInt8, RNA_B)   == UInt8(13)
                @test convert(UInt8, RNA_N)   == UInt8(14)
                @test convert(UInt8, RNA_Gap) == UInt8(15)
            end
        end

        @testset "UInt64" begin
            @testset "DNA conversions from UInt64" begin
                @test convert(DNANucleotide, UInt64(0)) == DNA_A
                @test convert(DNANucleotide, UInt64(1)) == DNA_C
                @test convert(DNANucleotide, UInt64(2)) == DNA_G
                @test convert(DNANucleotide, UInt64(3)) == DNA_T
                @test convert(DNANucleotide, UInt64(14)) == DNA_N
            end

            @testset "RNA conversions from UInt64" begin
                @test convert(RNANucleotide, UInt64(0)) == RNA_A
                @test convert(RNANucleotide, UInt64(1)) == RNA_C
                @test convert(RNANucleotide, UInt64(2)) == RNA_G
                @test convert(RNANucleotide, UInt64(3)) == RNA_U
                @test convert(RNANucleotide, UInt64(14)) == RNA_N
            end

            @testset "DNA conversions to UInt64" begin
                @test convert(UInt64, DNA_A) == UInt64(0)
                @test convert(UInt64, DNA_C) == UInt64(1)
                @test convert(UInt64, DNA_G) == UInt64(2)
                @test convert(UInt64, DNA_T) == UInt64(3)
                @test convert(UInt64, DNA_N) == UInt64(14)
            end

            @testset "RNA conversions to UInt64" begin
                @test convert(UInt64, RNA_A) == UInt64(0)
                @test convert(UInt64, RNA_C) == UInt64(1)
                @test convert(UInt64, RNA_G) == UInt64(2)
                @test convert(UInt64, RNA_U) == UInt64(3)
                @test convert(UInt64, RNA_N) == UInt64(14)
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
            @test convert(Int, DNA_N) == 14
            @test convert(DNANucleotide, 0) == DNA_A
            @test convert(DNANucleotide, 1) == DNA_C
            @test convert(DNANucleotide, 2) == DNA_G
            @test convert(DNANucleotide, 3) == DNA_T
            @test convert(DNANucleotide, 14) == DNA_N

            @test convert(Int, RNA_A) == 0
            @test convert(Int, RNA_C) == 1
            @test convert(Int, RNA_G) == 2
            @test convert(Int, RNA_U) == 3
            @test convert(Int, RNA_N) == 14
            @test convert(RNANucleotide, 0) == RNA_A
            @test convert(RNANucleotide, 1) == RNA_C
            @test convert(RNANucleotide, 2) == RNA_G
            @test convert(RNANucleotide, 3) == RNA_U
            @test convert(RNANucleotide, 14) == RNA_N
        end
    end

    @testset "Encoder" begin
        @testset "DNA" begin
            encode = Bio.Seq.encode
            EncodeError = Bio.Seq.EncodeError

            # 2 bits
            @test encode(DNAAlphabet{2}, DNA_A) === convert(UInt8, DNA_A) === 0b00
            @test encode(DNAAlphabet{2}, DNA_C) === convert(UInt8, DNA_C) === 0b01
            @test encode(DNAAlphabet{2}, DNA_G) === convert(UInt8, DNA_G) === 0b10
            @test encode(DNAAlphabet{2}, DNA_T) === convert(UInt8, DNA_T) === 0b11
            @test_throws EncodeError encode(DNAAlphabet{2}, DNA_M)
            @test_throws EncodeError encode(DNAAlphabet{2}, DNA_N)

            # 4 bits
            for nt in alphabet(DNANucleotide)
                @test encode(DNAAlphabet{4}, nt) === convert(UInt8, nt)
            end
            @test_throws EncodeError encode(DNAAlphabet{4}, Bio.Seq.DNA_INVALID)
        end
        @testset "RNA" begin
            encode = Bio.Seq.encode
            EncodeError = Bio.Seq.EncodeError

            # 2 bits
            @test encode(RNAAlphabet{2}, RNA_A) === convert(UInt8, RNA_A) === 0b00
            @test encode(RNAAlphabet{2}, RNA_C) === convert(UInt8, RNA_C) === 0b01
            @test encode(RNAAlphabet{2}, RNA_G) === convert(UInt8, RNA_G) === 0b10
            @test encode(RNAAlphabet{2}, RNA_U) === convert(UInt8, RNA_U) === 0b11
            @test_throws EncodeError encode(RNAAlphabet{2}, RNA_M)
            @test_throws EncodeError encode(RNAAlphabet{2}, RNA_N)

            # 4 bits
            for nt in alphabet(RNANucleotide)
                @test encode(RNAAlphabet{4}, nt) === convert(UInt8, nt)
            end
            @test_throws EncodeError encode(RNAAlphabet{4}, Bio.Seq.RNA_INVALID)
        end
    end

    @testset "Decoder" begin
        @testset "DNA" begin
            decode = Bio.Seq.decode
            DecodeError = Bio.Seq.DecodeError

            # 2 bits
            @test decode(DNAAlphabet{2}, 0b00) === convert(DNANucleotide, 0b00) === DNA_A
            @test decode(DNAAlphabet{2}, 0b01) === convert(DNANucleotide, 0b01) === DNA_C
            @test decode(DNAAlphabet{2}, 0b10) === convert(DNANucleotide, 0b10) === DNA_G
            @test decode(DNAAlphabet{2}, 0b11) === convert(DNANucleotide, 0b11) === DNA_T
            @test_throws DecodeError decode(DNAAlphabet{2}, 0b0100)
            @test_throws DecodeError decode(DNAAlphabet{2}, 0b1110)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(DNAAlphabet{4}, x) === convert(DNANucleotide, x)
            end
            @test_throws DecodeError decode(DNAAlphabet{4}, 0b10000)
        end
        @testset "RNA" begin
            decode = Bio.Seq.decode
            DecodeError = Bio.Seq.DecodeError

            # 2 bits
            @test decode(RNAAlphabet{2}, 0b00) === convert(RNANucleotide, 0b00) === RNA_A
            @test decode(RNAAlphabet{2}, 0b01) === convert(RNANucleotide, 0b01) === RNA_C
            @test decode(RNAAlphabet{2}, 0b10) === convert(RNANucleotide, 0b10) === RNA_G
            @test decode(RNAAlphabet{2}, 0b11) === convert(RNANucleotide, 0b11) === RNA_U
            @test_throws DecodeError decode(RNAAlphabet{2}, 0b0100)
            @test_throws DecodeError decode(RNAAlphabet{2}, 0b1110)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(RNAAlphabet{4}, x) === convert(RNANucleotide, x)
            end
            @test_throws DecodeError decode(RNAAlphabet{4}, 0b10000)
        end
    end

    @testset "Arithmetic and Order" begin
        @testset "DNA" begin
            @test DNA_A + 1 == DNA_C
            @test DNA_A + 2 == DNA_G
            @test DNA_T - 1 == DNA_G
            @test DNA_T - 2 == DNA_C
            @test DNA_T - DNA_A == 3
            @test DNA_T - DNA_C == 2
            @test DNA_A < DNA_C < DNA_G < DNA_T < DNA_M < DNA_N < DNA_Gap
            @test !(DNA_A > DNA_G)

            @test gap(DNANucleotide) === DNA_Gap
            @test collect(alphabet(DNANucleotide)) == [
                DNA_A, DNA_C, DNA_G, DNA_T,
                DNA_M, DNA_R, DNA_W, DNA_S,
                DNA_Y, DNA_K, DNA_V, DNA_H,
                DNA_D, DNA_B, DNA_N, DNA_Gap
            ]
        end
        @testset "RNA" begin
            @test RNA_A + 1 == RNA_C
            @test RNA_A + 2 == RNA_G
            @test RNA_U - 1 == RNA_G
            @test RNA_U - 2 == RNA_C
            @test RNA_U - RNA_A == 3
            @test RNA_U - RNA_C == 2
            @test RNA_A < RNA_C < RNA_G < RNA_U < RNA_M < RNA_N < RNA_Gap
            @test !(RNA_A > RNA_G)

            @test gap(RNANucleotide) === RNA_Gap
            @test collect(alphabet(RNANucleotide)) == [
                RNA_A, RNA_C, RNA_G, RNA_U,
                RNA_M, RNA_R, RNA_W, RNA_S,
                RNA_Y, RNA_K, RNA_V, RNA_H,
                RNA_D, RNA_B, RNA_N, RNA_Gap
            ]
        end
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
end

@testset "Aminoacids" begin
    @testset "Arithmetic and Order" begin
        @test AA_A + 1 == AA_R
        @test AA_R + 1 == AA_N
        @test AA_A + 2 == AA_N
        @test AA_R - 1 == AA_A
        @test AA_N - 2 == AA_A
        @test AA_D - AA_A ==  3
        @test AA_A - AA_D == -3
        @test (AA_A < AA_R < AA_N < AA_V < AA_O < AA_U <
               AA_B < AA_J < AA_Z < AA_X < AA_Term < AA_Gap)
        @test !(AA_J < AA_B)

        @test gap(AminoAcid) === AA_Gap
        @test length(alphabet(AminoAcid)) == 28
        @test AA_A in alphabet(AminoAcid)
        @test AA_I in alphabet(AminoAcid)
        @test AA_U in alphabet(AminoAcid)
    end

    @testset "Encoder" begin
        encode = Bio.Seq.encode
        @test encode(AminoAcidAlphabet, AA_A) === 0x00
        for aa in alphabet(AminoAcid)
            @test encode(AminoAcidAlphabet, aa) === convert(UInt8, aa)
        end
        @test_throws Bio.Seq.EncodeError encode(AminoAcidAlphabet, Bio.Seq.AA_INVALID)
    end

    @testset "Decoder" begin
        decode = Bio.Seq.decode
        @test decode(AminoAcidAlphabet, 0x00) === AA_A
        for x in 0x00:0x1b
            @test decode(AminoAcidAlphabet, x) === convert(AminoAcid, x)
        end
        @test_throws Bio.Seq.DecodeError decode(AminoAcidAlphabet, 0x1c)
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
                ("O", "PYL", AA_O),
                ("U", "SEC", AA_U),
                ("B", "ASX", AA_B),
                ("J", "XLE", AA_J),
                ("Z", "GLX", AA_Z),
                ("X", "XAA", AA_X),
            ]
            @test length(aas) == 26
            for (one, three, aa) in aas
                @test parse(AminoAcid, one) == aa
                @test parse(AminoAcid, three) == aa
            end
            @test parse(AminoAcid, "*") == AA_Term
            @test parse(AminoAcid, "-") == AA_Gap
        end

        @testset "Invalid Cases" begin
            @test_throws Exception parse(AminoAcid, "")
            @test_throws Exception parse(AminoAcid, "AL")
            @test_throws Exception parse(AminoAcid, "LA")
            @test_throws Exception parse(AminoAcid, "ALAA")
        end
    end
end

@testset "BioSequence" begin
    @testset "Constructing empty sequences" begin
        @test DNASequence() == BioSequence(DNANucleotide)
        @test RNASequence() == BioSequence(RNANucleotide)
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

        for len in [0, 1, 2, 3, 10, 32, 1000, 10000]
            test_string_construction(DNAAlphabet{4}, random_dna(len))
            test_string_construction(DNAAlphabet{4}, lowercase(random_dna(len)))
            test_string_construction(RNAAlphabet{4}, lowercase(random_rna(len)))
            test_string_construction(RNAAlphabet{4}, random_rna(len))
            test_string_construction(AminoAcidAlphabet, random_aa(len))
            test_string_construction(AminoAcidAlphabet, lowercase(random_aa(len)))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_string_construction(DNAAlphabet{2}, random_dna(len, probs))
            test_string_construction(DNAAlphabet{2}, lowercase(random_dna(len, probs)))
            test_string_construction(RNAAlphabet{2}, random_rna(len, probs))
            test_string_construction(RNAAlphabet{2}, lowercase(random_rna(len, probs)))
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
        EncodeError = Bio.Seq.EncodeError
        @test_throws EncodeError convert(BioSequence{DNAAlphabet{2}}, dna"AN")
        @test_throws EncodeError convert(BioSequence{RNAAlphabet{2}}, rna"AN")
    end

    @testset "Conversion between RNA and DNA" begin
        @test convert(RNASequence, DNASequence("ACGTN")) == rna"ACGUN"
        @test convert(DNASequence, RNASequence("ACGUN")) == dna"ACGTN"
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
            @test convert(ASCIIString, seq) == uppercase(str)
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
            @test convert(ASCIIString, seq) == uppercase(str)
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
                @test convert(ASCIIString, bioseq[part]) == seq[part]
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
    end

    @testset "Transformations" begin
        function test_reverse(A, seq)
            revseq = reverse(BioSequence{A}(seq))
            @test convert(ASCIIString, revseq) == reverse(seq)
        end

        function test_dna_complement(A, seq)
            comp = Seq.complement(BioSequence{A}(seq))
            @test convert(ASCIIString, comp) == dna_complement(seq)
        end

        function test_rna_complement(A, seq)
            comp = Seq.complement(BioSequence{A}(seq))
            @test convert(ASCIIString, comp) == rna_complement(seq)
        end

        function test_dna_revcomp(A, seq)
            revcomp = reverse_complement(BioSequence{A}(seq))
            @test convert(ASCIIString, revcomp) == reverse(dna_complement(seq))
        end

        function test_rna_revcomp(A, seq)
            revcomp = reverse_complement(BioSequence{A}(seq))
            @test convert(ASCIIString, revcomp) == reverse(rna_complement(seq))
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

    @testset "Mismatches" begin
        @test mismatches(dna"ACGT", dna"ACGT") == 0
        @test mismatches(dna"ACGT", dna"ACGTT") == 0
        @test mismatches(dna"ACGT", dna"ACGA") == 1
        @test mismatches(dna"ACGT", dna"ACGAA") == 1

        @test mismatches(rna"ACGU", rna"ACGU") == 0
        @test mismatches(rna"ACGU", rna"ACGUU") == 0
        @test mismatches(rna"ACGU", rna"ACGA") == 1
        @test mismatches(rna"ACGU", rna"ACGAA") == 1

        @test mismatches(aa"MTTQAP", aa"MTTQAP") == 0
        @test mismatches(aa"MTTQAP", aa"MTTQAPM") == 0
        @test mismatches(aa"MTTQAP", aa"MTTQAT") == 1
        @test mismatches(aa"MTTQAP", aa"MTTQATT") == 1

        @test mismatches(dna"ACGT", dna"TACTG"[2:end]) == 2
        @test mismatches(dna"ACGT"[2:end], dna"AGT") == 1

        function test_mismatches(A, a, b)
            count = 0
            for (ca, cb) in zip(a, b)
                if ca != cb
                    count += 1
                end
            end
            seq_a = BioSequence{A}(a)
            seq_b = BioSequence{A}(b)
            @test mismatches(seq_a, seq_b) === mismatches(seq_b, seq_a) == count
        end

        for len in [0, 1, 10, 32, 1000, 10000, 100000], _ in 1:10
            test_mismatches(DNAAlphabet{4}, random_dna(len), random_dna(len))
            test_mismatches(RNAAlphabet{4}, random_rna(len), random_rna(len))
            test_mismatches(AminoAcidAlphabet, random_aa(len), random_aa(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_mismatches(DNAAlphabet{2}, random_dna(len, probs), random_dna(len, probs))
            test_mismatches(RNAAlphabet{2}, random_rna(len, probs), random_rna(len, probs))
        end
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
end

@testset "Composition" begin
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
        seq_counts = composition(DNASequence(seq))
        return string_counts[DNA_A] == seq_counts[DNA_A] &&
               string_counts[DNA_C] == seq_counts[DNA_C] &&
               string_counts[DNA_G] == seq_counts[DNA_G] &&
               string_counts[DNA_T] == seq_counts[DNA_T] &&
               string_counts[DNA_N] == seq_counts[DNA_N]
    end

    function check_nucleotide_count(::Type{RNANucleotide}, seq::AbstractString)
        string_counts = string_nucleotide_count(RNANucleotide, seq)
        seq_counts = composition(RNASequence(seq))
        return string_counts[RNA_A] == seq_counts[RNA_A] &&
               string_counts[RNA_C] == seq_counts[RNA_C] &&
               string_counts[RNA_G] == seq_counts[RNA_G] &&
               string_counts[RNA_U] == seq_counts[RNA_U] &&
               string_counts[RNA_N] == seq_counts[RNA_N]
    end

    function check_kmer_nucleotide_count(::Type{DNANucleotide}, seq::AbstractString)
        string_counts = string_nucleotide_count(DNANucleotide, seq)
        kmer_counts = composition(dnakmer(seq))
        return string_counts[DNA_A] == kmer_counts[DNA_A] &&
               string_counts[DNA_C] == kmer_counts[DNA_C] &&
               string_counts[DNA_G] == kmer_counts[DNA_G] &&
               string_counts[DNA_T] == kmer_counts[DNA_T] &&
               string_counts[DNA_N] == kmer_counts[DNA_N]
    end

    function check_kmer_nucleotide_count(::Type{RNANucleotide}, seq::AbstractString)
        string_counts = string_nucleotide_count(RNANucleotide, seq)
        kmer_counts = composition(rnakmer(seq))
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

    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_E] === 0
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_A] === 1
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_P] === 2
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_V] === 3

    for len in [1, 10, 32]
        @test all(Bool[check_kmer_nucleotide_count(DNANucleotide, random_dna_kmer(len)) for _ in 1:reps])
        @test all(Bool[check_kmer_nucleotide_count(RNANucleotide, random_rna_kmer(len)) for _ in 1:reps])
    end
end

@testset "Kmer" begin
    reps = 10
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

        # Check that kmers can be constructed from a BioSequence
        #   BioSequence → Kmer → BioSequence
        function check_biosequence_construction(seq::BioSequence)
            return convert(BioSequence, convert(Kmer, seq)) == seq
        end

        # Check that kmers can be constructed from an array of nucleotides
        #   Vector{Nucleotide} → Kmer → Vector{Nucleotide}
        function check_nucarray_kmer{T <: Nucleotide}(seq::Vector{T})
            return convert(AbstractString, [convert(Char, c) for c in seq]) ==
                   convert(AbstractString, kmer(seq...))
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
            @test all(Bool[check_uint64_convertion(DNANucleotide, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])
            @test all(Bool[check_uint64_convertion(RNANucleotide, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])

            # String construction
            @test all(Bool[check_string_construction(DNANucleotide, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(RNANucleotide, random_rna_kmer(len)) for _ in 1:reps])

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
            @test dnakmer("ACG") != rnakmer("ACG")
            @test dnakmer("T")   != rnakmer("U")
            @test dnakmer("AC")  != dnakmer("AG")
            @test rnakmer("AC")  != rnakmer("AG")

            @test dnakmer("ACG") != rna"ACG"
            @test dnakmer("T")   != rna"U"
            @test dnakmer("AC")  != dna"AG"
            @test rnakmer("AC")  != rna"AG"

            @test rna"ACG" != dnakmer("ACG")
            @test rna"U"   != dnakmer("T")
            @test dna"AG"  != dnakmer("AC")
            @test rna"AG"  != rnakmer("AC")
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
            @test start(dnakmer("ACTG")) == 1
            @test start(dnakmer(""))     == 1

            @test next(dnakmer("ACTG"), 1) == (DNA_A, 2)
            @test next(dnakmer("ACTG"), 4) == (DNA_G, 5)

            @test  done(dnakmer(""), 1)
            @test !done(dnakmer("ACTG"), 1)
            @test !done(dnakmer("ACTG"), 4)
            @test  done(dnakmer("ACTG"), 5)
            @test !done(dnakmer("ACTG"), -1)


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
            @test start(rnakmer("ACUG")) == 1
            @test start(rnakmer(""))     == 1

            @test next(rnakmer("ACUG"), 1) == (RNA_A, 2)
            @test next(rnakmer("ACUG"), 4) == (RNA_G, 5)

            @test  done(rnakmer(""), 1)
            @test !done(rnakmer("ACUG"), 1)
            @test !done(rnakmer("ACUG"), 4)
            @test  done(rnakmer("ACUG"), 5)
            @test !done(rnakmer("ACUG"), -1)

            rna_vec = [RNA_A, RNA_C, RNA_U, RNA_G]
            @test all([nt === rna_vec[i] for (i, nt) in enumerate(rna_kmer)])
        end
    end

    @testset "Find" begin
        kmer = dnakmer("ACGAG")

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

    @testset "Transformations" begin
        function test_reverse(T, seq)
            revseq = reverse(Kmer{T,length(seq)}(seq))
            @test convert(ASCIIString, revseq) == reverse(seq)
        end

        function test_dna_complement(seq)
            comp = Seq.complement(DNAKmer{length(seq)}(seq))
            @test convert(ASCIIString, comp) == dna_complement(seq)
        end

        function test_rna_complement(seq)
            comp = Seq.complement(RNAKmer{length(seq)}(seq))
            @test convert(ASCIIString, comp) == rna_complement(seq)
        end

        function test_dna_revcomp(seq)
            revcomp = reverse_complement(DNAKmer{length(seq)}(seq))
            @test convert(ASCIIString, revcomp) == reverse(dna_complement(seq))
        end

        function test_rna_revcomp(seq)
            revcomp = reverse_complement(RNAKmer{length(seq)}(seq))
            @test convert(ASCIIString, revcomp) == reverse(rna_complement(seq))
        end

        @testset "Reverse" begin
            for len in 0:32, _ in 1:10
                test_reverse(DNANucleotide, random_dna_kmer(len))
                test_reverse(RNANucleotide, random_rna_kmer(len))
            end
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
            test_mismatches(dnakmer(a), dnakmer(b))

            a = random_rna_kmer(len)
            b = random_rna_kmer(len)
            test_mismatches(rnakmer(a), rnakmer(b))
        end
    end

    @testset "EachKmer" begin
        function string_eachkmer(seq::AbstractString, k, step)
            kmers = ASCIIString[]
            i = 1
            for i in 1:step:length(seq) - k + 1
                subseq = seq[i:i + k - 1]
                if !in('N', subseq)
                    push!(kmers, subseq)
                end
            end
            return kmers
        end

        function test_eachkmer(A, seq::AbstractString, k, step)
            T = eltype(A)
            xs = [convert(AbstractString, x) for (i, x) in collect(each(Kmer{T,k}, BioSequence{A}(seq), step))]
            ys = [convert(AbstractString, x) for (i, x) in collect(eachkmer(BioSequence{A}(seq), k, step))]
            zs = string_eachkmer(seq, k, step)
            @test xs == ys == zs
        end

        len = 10000

        for k in [0, 1, 3, 16, 32], step in 1:3, _ in 1:10
            test_eachkmer(DNAAlphabet{4}, random_dna(len), k, step)
            test_eachkmer(RNAAlphabet{4}, random_rna(len), k, step)

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_eachkmer(DNAAlphabet{2}, random_dna(len, probs), k, step)
            test_eachkmer(RNAAlphabet{2}, random_rna(len, probs), k, step)
        end

        @test isempty(collect(each(DNAKmer{1}, dna"")))
        @test isempty(collect(each(DNAKmer{1}, dna"NNNNNNNNNN")))
        @test_throws Exception each(DNAKmer{-1}, dna"ACGT")
        @test_throws Exception each(DNAKmer{33}, dna"ACGT")
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

        function test_kmer_count(A, seq::AbstractString, k, step)
            T = eltype(A)
            string_counts = string_kmer_count(T, seq, k, step)
            kmer_counts = KmerCounts{T,k}(BioSequence{A}(seq), step)
            ok = true
            for y in UInt64(0):UInt64(4^k-1)
                x = convert(Kmer{T,k}, y)
                if string_counts[x] != kmer_counts[x]
                    ok = false
                    break
                end
            end
            @test ok
        end

        for len in [1, 10, 32, 1000, 10000], k in [1, 2, 5], step in 1:3, _ in 1:10
            test_kmer_count(DNAAlphabet{4}, random_dna(len), k, step)
            test_kmer_count(RNAAlphabet{4}, random_rna(len), k, step)

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_kmer_count(DNAAlphabet{2}, random_dna(len, probs), k, step)
            test_kmer_count(RNAAlphabet{2}, random_rna(len, probs), k, step)
        end
    end
end

@testset "Translation" begin
    # crummy string translation to test against
    standard_genetic_code_dict = Dict{ASCIIString,Char}(
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


@testset "Parsing" begin
    @testset "FASTA" begin
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
            if contains(tags, "comments")
                # currently comments are not supported
                continue
            end
            if valid
                @test check_fasta_parse(joinpath(path, specimen["filename"]))
            else
                @test_throws Exception check_fasta_parse(joinpath(path, specimen["filename"]))
            end
        end
    end

    @testset "FASTQ" begin
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
end

end # TestSeq
