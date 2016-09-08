module TestVar

using Base.Test

using Bio: Seq, Var
using TestFunctions


# A naieve site counting function which will also compute the answers to make
# sure the clever bit parallel methods are giving the right answers...

@inline function iscase{T<:Nucleotide}(::Type{Ambiguous}, a::T)
    return isambiguous(a)
end

@inline function iscase{T<:Nucleotide}(::Type{Ambiguous}, a::T, b::T)
    return isambiguous(a) | isambiguous(b)
end

@inline function iscase{T<:Nucleotide}(::Type{Gap}, a::T)
    return reinterpret(UInt8, a) == 0
end

@inline function iscase{T<:Nucleotide}(::Type{Gap}, a::T, b::T)
    return (reinterpret(UInt8, a) == 0) | (reinterpret(UInt8, b) == 0)
end

@inline function iscase{T<:Nucleotide}(::Type{Pairdel}, a::T, b::T)
    return iscase(Ambiguous, a, b) | iscase(Gap, a, b)
end

@inline function iscase{T<:Nucleotide}(::Type{Mutated}, a::T, b::T)
    return !iscase(Pairdel, a, b) & (a != b)
end

@inline function iscase{T<:Nucleotide}(::Type{Conserved}, a::T, b::T)
    return !iscase(Pairdel, a, b) & (a == b)
end

@inline function iscase{T<:Nucleotide}(::Type{Transition}, a::T, b::T)
    return iscase(Mutated, a, b) & ((ispurine(a) & ispurine(b)) | (ispyrimidine(a) & ispyrimidine(b)))
end

@inline function iscase{T<:Nucleotide}(::Type{Transversion}, a::T, b::T)
    return iscase(Mutated, a, b) & ((ispurine(a) & ispyrimidine(b)) | (ispyrimidine(a) & ispurine(b)))
end

function count_sites_naieve{T<:SiteCase,A<:Alphabet}(::Type{T}, a::BioSequence{A}, b::BioSequence{A})
    if length(a) != length(b)
        error("`a` and `b` must be the same length.")
    end
    n = 0
    @inbounds for (i, j) in zip(a, b)
        n += ifelse(iscase(T, i, j), 1, 0)
    end
    return n
end

function count_sites_naieve{T<:Union{Gap,Ambiguous}}(::Type{T}, a::BioSequence{DNAAlphabet{4}})
    n = 0
    @inbounds for i in a
        n += ifelse(iscase(T, i), 1, 0)
    end
    return n
end

function generate_testcase{A<:Union{DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4}, RNAAlphabet{2}}}(::Type{A}, len::Int)
    a = [convert(Char, i)  for i in alphabet(A)]
    probs = Vector{Float64}(length(a))
    fill!(probs, 1 / length(a))
    return BioSequence{A}(random_seq(len, a, probs))
end


@testset "Var" begin

    @testset "Site counting and identification" begin

        @testset "Internals" begin

            @testset "Bit parallel operations" begin

                # 4 consecutive Conservedes, 4 consecutive misConservedes, 4 consecutive ambigs,
                # 4 consecutive gaps.
                aseq = dna"ATCGATCGARMGA-M-"
                bseq = dna"ATCGTCGARMWY-T--"
                cseq = dna"ATCGMRWSYKVHDBN-"
                a = aseq.data[1]
                b = bseq.data[1]
                c = cseq.data[1]

                @testset "Counting zeros" begin
                    @test Var.count_zero_nibbles(0x0000000000000000) == 16
                    @test Var.count_zero_nibbles(0xF004020000403010) == 10
                    @test Var.count_zero_nibbles(0xFFFFFFFFFFFFFFFF) == 0
                end

                @testset "Counting ones" begin
                    @test Var.count_one_nibbles(0x01011C1111F11011) == 1
                    @test Var.count_one_nibbles(0xF4011C1111F1101F) == 3
                    @test Var.count_one_nibbles(0x0000000000000000) == 0
                    @test Var.count_one_nibbles(0xFFFFFFFFFFFFFFFF) == 16
                end

                @testset "Enumerating nibbles" begin
                    @test Var.enumerate_nibbles(0x0000000000000000) == 0x0000000000000000
                    @test Var.enumerate_nibbles(0xF004020000403010) == 0x4001010000102010
                end

                @testset "Masking nibbles" begin
                    @test Var.create_nibble_mask(Gap, a) == 0xF0F0000000000000
                    @test Var.create_nibble_mask(Gap, b) == 0xFF0F000000000000
                    @test Var.create_nibble_mask(Gap, a, b) == 0xFFFF000000000000
                    @test Var.create_nibble_mask(Ambiguous, a) == 0x0F000FF000000000
                    @test Var.create_nibble_mask(Ambiguous, b) == 0x0000FFFF00000000
                    @test Var.create_nibble_mask(Ambiguous, a, b) == 0x0F00FFFF00000000
                    @test Var.create_nibble_mask(Pairdel, a) == 0xFFF00FF000000000
                    @test Var.create_nibble_mask(Pairdel, b) == 0xFF0FFFFF00000000
                    @test Var.create_nibble_mask(Pairdel, a, b) == 0xFFFFFFFF00000000
                    @test Var.create_nibble_mask(Conserved, a, b) == 0x000000000000FFFF
                    @test Var.create_nibble_mask(Mutated, a, b) == 0x00000000FFFF0000
                end

                @testset "Case counting" begin

                    @testset "Gaps" begin
                        # Set cases we've reasoned about to test certain
                        # behaviours and properties.
                        @test Var.count_sites4(Gap, a) == 2
                        @test Var.count_sites4(Gap, b) == 3
                        @test Var.count_sites4(Gap, a | b) == Var.count_sites4(Gap, b | a) == 1
                        @test Var.count_sites4(Gap, a, b) == Var.count_sites4(Gap, b, a) == 4
                        @test Var.count_sites4(Gap, a | c) == Var.count_sites4(Gap, c | a) == 1
                        @test Var.count_sites4(Gap, a, c) == Var.count_sites4(Gap, c, a) == 2
                        @test Var.count_sites4(Gap, b, c) == Var.count_sites4(Gap, c, b) == 3
                        @test Var.count_sites4(Gap, b | c) == Var.count_sites4(Gap, c | b) == 1
                        # Randomly generated cases.
                        for i in 1:100
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Gap, s)
                            off = 16 - n
                            @test (Var.count_sites4(Gap, s.data[1]) - off) == expected
                        end
                    end

                    @testset "Ambiguities" begin
                        # Set cases we've reasoned about to test certain
                        # behaviours and properties.
                        @test Var.count_sites4(Ambiguous, c) == 11
                        @test Var.count_sites4(Ambiguous, a) == 3
                        @test Var.count_sites4(Ambiguous, b) == 4
                        @test Var.count_sites4(Ambiguous, a, b) == Var.count_sites4(Ambiguous, b, a) == 5
                        @test Var.count_sites4(Ambiguous, a, c) == Var.count_sites4(Ambiguous, c, a) == 11
                        @test Var.count_sites4(Ambiguous, b, c) == Var.count_sites4(Ambiguous, c, b) == 11
                        # Randomly generated cases.
                        for i in 1:100
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected_a = count_sites_naieve(Ambiguous, s)
                            expected_b = count_sites_naieve(Ambiguous, s, s2)
                            @test Var.count_sites4(Ambiguous, s.data[1]) == expected_a
                            @test Var.count_sites4(Ambiguous, s.data[1], s2.data[1]) == expected_b
                        end
                    end

                    @testset "Pairdel" begin
                        # Set cases we've reasoned about to test certain
                        # behaviours and properties.
                        @test Var.count_sites4(Pairdel, a) == 5
                        @test Var.count_sites4(Pairdel, b) == 7
                        @test Var.count_sites4(Pairdel, c) == 12
                        @test Var.count_sites4(Pairdel, a, b) == Var.count_sites4(Pairdel, b, a) == 8
                        @test Var.count_sites4(Pairdel, a, c) == Var.count_sites4(Pairdel, c, a) == 12
                        @test Var.count_sites4(Pairdel, b, c) == Var.count_sites4(Pairdel, c, b) == 12
                        # Randomly generated testcases.
                        for i in 1:100
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Pairdel, s, s2)
                            off = 16 - n
                            @test (Var.count_sites4(Pairdel, s.data[1], s2.data[1]) - off) == expected
                        end
                    end

                    @testset "Conserved" begin
                        # Set cases we've reasoned about to test certain
                        # behaviours and properties.
                        @test Var.count_sites4(Conserved, a, b) == Var.count_sites4(Conserved, b, a) == 4
                        @test Var.count_sites4(Conserved, a, c) == Var.count_sites4(Conserved, c, a) == 4
                        @test Var.count_sites4(Conserved, b, c) == Var.count_sites4(Conserved, c, b) == 4
                        # Randomly generated cases.
                        for i in 1:100
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Conserved, s, s2)
                            @test Var.count_sites4(Conserved, s.data[1], s2.data[1]) == expected
                        end
                    end

                    @testset "Mutated" begin
                        # Set cases we've reasoned about to test certain
                        # behaviours and properties.
                        @test Var.count_sites4(Mutated, a, b) == Var.count_sites4(Mutated, b, a) == 4
                        @test Var.count_sites4(Mutated, a, c) == Var.count_sites4(Mutated, c, a) == 0
                        @test Var.count_sites4(Mutated, b, c) == Var.count_sites4(Mutated, c, b) == 0
                        # Randomly generated cases.
                        for i in 1:100
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Mutated, s, s2)
                            @test Var.count_sites4(Mutated, s.data[1], s2.data[1]) == expected
                        end
                    end

                    @testset "Transition" begin
                        # Set cases we've reasoned about to test certain
                        # behaviours and properties.
                        @test Var.count_sites4(Transition, a, b) == Var.count_sites4(Transition, b, a) == 2
                        @test Var.count_sites4(Transition, a, c) == Var.count_sites4(Transition, b, c) == 0
                        @test Var.count_sites4(Transition, b, c) == Var.count_sites4(Transition, c, b) == 0
                        # Randomly generated cases.
                        for i in 1:100
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Transition, s, s2)
                            @test Var.count_sites4(Transition, s.data[1], s2.data[1]) == expected
                        end
                    end
                end
            end

            @testset "ShiftedIntsItr" begin
                @testset "4 bit biological sequences" begin

                    # Replace test cases with automatic generation of any shit.


                    tinyseq = dna"ATCG"
                    smallseq = dna"AAAATTTTGGGGCCCC"
                    largeseq = dna"AAAATTTTGGGGCCCCAAAATTTTGGGGCCCC"
                    largerseq = dna"AAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAAGGGG"

                    # Very basic cases - no shifting.
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq)) == tinyseq.data
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(smallseq)) == smallseq.data
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq)) == largeseq.data
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largerseq)) == largerseq.data

                    # Subsetting cases - tinyseq
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq[1:4])) == tinyseq.data
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq[2:4])) == [0x0000000000000428]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq[3:4])) == [0x0000000000000042]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq[4:4])) == [0x0000000000000004]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq[1:3])) == [0x0000000000000281]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq[2:3])) == [0x0000000000000028]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(tinyseq[3:3])) == [0x0000000000000002]

                    # Subsetting cases - smallseq
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(smallseq[1:16])) == smallseq.data
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(smallseq[1:10])) == [0x0000004488881111]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(smallseq[3:10])) == [0x0000000044888811]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(smallseq[5:16])) == [0x0000222244448888]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(smallseq[5:10])) == [0x0000000000448888]

                    # Subsetting cases - largeseq
                    ## Subsetting cases that span one integer.
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq)) == largeseq.data
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[1:16])) == [0x2222444488881111]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[1:5])) == [0x0000000000081111]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[5:10])) == [0x0000000000448888]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[5:16])) == [0x0000222244448888]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[17:32])) == [0x2222444488881111]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[17:28])) == [0x0000444488881111]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[22:32])) == [0x0000022224444888]

                    ## Subsetting cases that span two integers.
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[5:32])) == [0x1111222244448888,
                                                                                       0x0000222244448888]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[5:20])) == [0x1111222244448888]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[5:21])) == [0x1111222244448888,
                                                                                       0x0000000000000008]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[6:21])) == [0x8111122224444888]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[11:21])) == [0x0000081111222244]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largeseq[17:21])) == [0x0000000000081111]

                    # Subsetting cases - largerseq
                    # Subsetting cases that span three integers.
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largerseq[5:40])) == [0x1111222244448888,
                                                                                        0x1111222244448888,
                                                                                        0x0000000000004444]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largerseq[8:39])) == [0x8881111222244448,
                                                                                        0x4441111222244448]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largerseq[12:30])) == [0x4448888111122224,
                                                                                         0x0000000000000224]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largerseq[3:37])) == [0x1122224444888811,
                                                                                        0x1122224444888811,
                                                                                        0x0000000000000411]
                    @test collect(Var.ShiftedInts{DNAAlphabet{4}}(largerseq[3:35])) == [0x1122224444888811,
                                                                                        0x1122224444888811,
                                                                                        0x0000000000000001]
                end
            end

            @testset "Looped case counting" begin
                @testset "Aligned sequence data" begin
                    # Create a 20bp test DNA sequence pair containing every possible transition (4),
                    # every possible transversion (8), and 2 gapped sites and 2 ambiguous sites.
                    # This leaves 4 sites non-mutated/conserved.

                    # Repeating that sequence pair 100 times gives a 2000 base
                    # sequence pair. With 400 tranitions, 800 transversions,
                    # 200 gapped sites and 200 ambiguous sites.
                    # 1200 Mutations, 400 conserved sites, 400 pairdel sites.
                    seq1 = dna"ATTG-ACCTGGNTTTCCGAA"^100
                    seq2 = dna"A-ACAGAGTATACRGTCGTC"^100

                    @testset "Verifying naieve methods" begin
                        @test count_sites_naieve(Gap, seq1, seq2) == 200
                        @test count_sites_naieve(Ambiguous, seq1, seq2) == 200
                        @test count_sites_naieve(Pairdel, seq1, seq2) == 400
                        @test count_sites_naieve(Conserved, seq1, seq2) == 400
                        @test count_sites_naieve(Mutated, seq1, seq2) == 1200
                        @test count_sites_naieve(Transition, seq1, seq2) == 400
                        @test count_sites_naieve(Transversion, seq1, seq2) == 800
                    end


                    #@test Var.count_sites4(Ambiguous, seq1.data, seq2.data) == 200
                    #@test Var.count_sites4(Conserved, seq1.data, seq2.data) == 400
                    #@test Var.count_sites4(Mutated, seq1.data, seq2.data) == 1200
                    #@test Var.count_sites4(Transition, seq1.data, seq2.data) == 400
                    #@test Var.count_sites4(Transversion, seq1.data, seq2.data) == 800

                    #@test Var.count_sites4(Ambiguous, seq1.data, seq2.data) == 200
                    #@test Var.count_sites4(Conserved, seq1.data, seq2.data) == 400
                    #@test Var.count_sites4(Mutated, seq1.data, seq2.data) == 1200
                    #@test Var.count_sites4(Transition, seq1.data, seq2.data) == 400
                end
            end
        end


#=
        @testset "Counting sites" begin
            @testset "One sequence" begin
                testseq1 = dna"AAAAATTTTTRM--YGGGGG"
                testseq2 = dna"ATCGNYR"
                testseq3 = dna"ATCGB"

                #@test count_sites(Ambiguous, testseq1) == 3
                #@test count_sites(Ambiguous, testseq2) == 3
                #@test count_sites(Ambiguous, testseq3) == 1
            end

            @testset "Two sequences" begin






            end

        end
=#
    end



#=
    @testset "Distance Computation" begin

        dnas1 = [dna"ATTG-ACCTGGNTTTCCGAA", dna"A-ACAGAGTATACRGTCGTC"]
        m1 = seqmatrix(dnas1, :seq)

        dnas2 = [dna"attgaacctggntttccgaa",
                 dna"atacagagtatacrgtcgtc"]
        dnas3 = [dna"attgaacctgtntttccgaa",
                 dna"atagaacgtatatrgccgtc"]
        m2 = seqmatrix(dnas2, :seq)

        @test distance(Count{AnyMutation}, dnas1) == ([12], [16])
        @test distance(Count{TransitionMutation}, dnas1) == ([4], [16])
        @test distance(Count{TransversionMutation}, dnas1) == ([8], [16])
        @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])
        @test distance(Count{AnyMutation}, m1) == ([12], [16])
        @test distance(Count{TransitionMutation}, m1) == ([4], [16])
        @test distance(Count{TransversionMutation}, m1) == ([8], [16])
        @test distance(Count{Kimura80}, m1) == ([4], [8], [16])

        @test distance(Count{AnyMutation}, dnas2, 5, 5)[1][:] == [2, 4, 3, 3]
        @test distance(Count{AnyMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
        @test distance(Count{TransitionMutation}, dnas2, 5, 5)[1][:] == [0, 2, 1, 1]
        @test distance(Count{TransitionMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
        @test distance(Count{TransversionMutation}, dnas2, 5, 5)[1][:] == [2, 2, 2, 2]
        @test distance(Count{TransversionMutation}, dnas2, 5, 5)[2][:] == [5, 5, 3, 5]
        @test distance(Count{Kimura80}, dnas1) == ([4], [8], [16])

        @test distance(Count{AnyMutation}, dnas2) == ([12], [18])
        @test distance(Count{TransitionMutation}, dnas2) == ([4], [18])
        @test distance(Count{TransversionMutation}, dnas2) == ([8], [18])
        @test distance(Count{Kimura80}, dnas2) == ([4], [8], [18])
        @test distance(Count{AnyMutation}, m2) == ([12], [18])
        @test distance(Count{TransitionMutation}, m2) == ([4], [18])
        @test distance(Count{TransversionMutation}, m2) == ([8], [18])
        @test distance(Count{Kimura80}, m2) == ([4], [8], [18])

        d = distance(Proportion{AnyMutation}, dnas2, 5, 5)
        a = [0.4, 0.8, 1.0, 0.6]
        for i in 1:length(d[1])
            @test_approx_eq_eps d[1][i] a[i] 1e-4
        end
        @test d[2][:] == [5, 5, 3, 5]
        d = distance(Proportion{TransitionMutation}, dnas2, 5, 5)
        a = [0.0, 0.4, 0.333333, 0.2]
        for i in 1:length(d[1])
            @test_approx_eq_eps d[1][i] a[i] 1e-4
        end
        @test d[2][:] == [5, 5, 3, 5]
        d = distance(Proportion{TransversionMutation}, dnas2, 5, 5)
        a = [0.4, 0.4, 0.666667, 0.4]
        for i in 1:length(d[1])
            @test_approx_eq_eps d[1][i] a[i] 1e-4
        end
        @test d[2][:] == [5, 5, 3, 5]

        @test distance(Proportion{AnyMutation}, dnas1) == ([(12 / 16)], [16])
        @test distance(Proportion{TransitionMutation}, dnas1) == ([(4 / 16)], [16])
        @test distance(Proportion{TransversionMutation}, dnas1) == ([(8 / 16)], [16])
        @test distance(Proportion{AnyMutation}, m1) == ([(12 / 16)], [16])
        @test distance(Proportion{TransitionMutation}, m1) == ([(4 / 16)], [16])
        @test distance(Proportion{TransversionMutation}, m1) == ([(8 / 16)], [16])

        @test distance(Proportion{AnyMutation}, dnas2) == ([(12 / 18)], [18])
        @test distance(Proportion{TransitionMutation}, dnas2) == ([(4 / 18)], [18])
        @test distance(Proportion{TransversionMutation}, dnas2) == ([(8 / 18)], [18])
        @test distance(Proportion{AnyMutation}, m2) == ([(12 / 18)], [18])
        @test distance(Proportion{TransitionMutation}, m2) == ([(4 / 18)], [18])
        @test distance(Proportion{TransversionMutation}, m2) == ([(8 / 18)], [18])

        @test distance(JukesCantor69, dnas1) == ([Inf], [Inf]) # Returns infinity as 12/16 is 0.75 - mutation saturation.
        @test distance(JukesCantor69, m1) == ([Inf], [Inf])

        @test round(distance(JukesCantor69, dnas2)[1][1], 3) == 1.648
        @test round(distance(JukesCantor69, dnas2)[2][1], 3) == 1
        @test round(distance(JukesCantor69, m2)[1][1], 3) == 1.648
        @test round(distance(JukesCantor69, m2)[2][1], 3) == 1
        @test_throws DomainError distance(JukesCantor69, dnas2, 5, 5)
        d = distance(JukesCantor69, dnas3, 5, 5)
        a = [0.232616, 0.571605, 0.44084, 0.571605]
        v = [0.0595041, 0.220408, 0.24, 0.220408]
        for i in 1:length(d[1])
            @test_approx_eq_eps d[1][i] a[i] 1e-5
            @test_approx_eq_eps d[2][i] v[i] 1e-5
        end

        @test round(distance(Kimura80, dnas2)[1][1], 3) == 1.648
        @test round(distance(Kimura80, dnas2)[2][1], 3) == 1
        @test round(distance(Kimura80, m2)[1][1], 3) == 1.648
        @test round(distance(Kimura80, m2)[2][1], 3) == 1

    end
    =#

end



end # module TestVar
