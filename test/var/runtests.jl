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

@inline function iscase{T<:Nucleotide}(::Type{Indel}, a::T)
    return reinterpret(UInt8, a) == 0
end

@inline function iscase{T<:Nucleotide}(::Type{Indel}, a::T, b::T)
    return (reinterpret(UInt8, a) == 0) | (reinterpret(UInt8, b) == 0)
end

@inline function iscase{T<:Nucleotide}(::Type{Certain}, a::T, b::T)
    return !iscase(Ambiguous, a, b) & !iscase(Indel, a, b)
end

@inline function iscase{T<:Nucleotide}(::Type{Match}, a::T, b::T)
    return a == b
end

@inline function iscase{T<:Nucleotide}(::Type{Mismatch}, a::T, b::T)
    return a != b
end

@inline function iscase{T<:Nucleotide}(::Type{Mutated}, a::T, b::T)
    return iscase(Certain, a, b) & (a != b)
end

@inline function iscase{T<:Nucleotide}(::Type{Conserved}, a::T, b::T)
    return iscase(Certain, a, b) & (a == b)
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

function count_sites_naieve{T<:Union{Indel,Ambiguous}}(::Type{T}, a::BioSequence{DNAAlphabet{4}})
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

function generate_masking_tester{A<:Union{DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4}, RNAAlphabet{2}}}(::Type{A})
    symbols = alphabet(A)
    arra = Vector{eltype(A)}()
    arrb = Vector{eltype(A)}()
    for i in 1:length(symbols), j in i:length(symbols)
        push!(arra, symbols[i])
        push!(arrb, symbols[j])
    end
    return BioSequence{A}(arra), BioSequence{A}(arrb)
end


@testset "Var" begin

    @testset "Site counting and identification" begin

        @testset "Internals" begin

            @testset "Nibble operations" begin

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

                    tester_a, tester_b = generate_masking_tester(DNAAlphabet{4})
                    a = tester_a.data
                    b = tester_b.data

                    gap_a_answers = [0xFFFFFFFFFFFFFFFF, 0x0000000000000000, 0xFFFFFFFF00000000]
                    gap_b_answers = [0x000000000000000F, 0x0000000000000000, 0xFFFFFFFF00000000]
                    gap_ab_answers = [0xFFFFFFFFFFFFFFFF, 0x0000000000000000, 0xFFFFFFFF00000000]

                    for i in 1:endof(a)
                        if i > 1
                            ans = i == endof(a) ? 3 : 2
                        else
                            ans = 1
                        end
                        @test Var.create_nibble_mask(Indel, a[i]) == gap_a_answers[ans]
                        @test Var.create_nibble_mask(Indel, b[i]) == gap_b_answers[ans]
                        @test Var.create_nibble_mask(Indel, a[i], b[i]) == gap_ab_answers[ans]
                    end

                    ambig_a_answers = [0x0000000000000000, 0x0000000000000000,
                                       0xFFF0000000000000, 0x000000FFFFFFFFFF,
                                       0xFFFFFFFFFF000000, 0xFFFFFFFFFFFFFFFF,
                                       0xFFFF00000000FFFF, 0xFFFFFFFFFFFFFFFF,
                                       0x00000000FFFFFFFF]

                    ambig_b_answers = [0xFFFFFFF0FFF0F000, 0x0FFFFFFF0FFF0F00,
                                       0xF0FFFFFFFF0FFF0F, 0xF0FFF0FFFFFFF0FF,
                                       0xFFFFFF0FFFFFFFFF, 0xFFF0FFFFFFFF0FFF,
                                       0xFFFFFFFFFFF0FFFF, 0xFFFFFFFFFFFFFFFF,
                                       0x00000000FFFFFFFF]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Ambiguous, a[i]) == ambig_a_answers[i]
                        @test Var.create_nibble_mask(Ambiguous, b[i]) == ambig_b_answers[i]
                        @test Var.create_nibble_mask(Ambiguous, a[i], b[i]) == (ambig_a_answers[i] | ambig_b_answers[i])
                    end

                    certain_a_answers = [0x0000000000000000, 0xFFFFFFFFFFFFFFFF,
                                         0x000FFFFFFFFFFFFF, 0xFFFFFF0000000000,
                                         0x0000000000FFFFFF, 0x0000000000000000,
                                         0x0000FFFFFFFF0000, 0x0000000000000000,
                                         0x0000000000000000]

                    certain_b_answers = [0x0000000F000F0FF0, 0xF0000000F000F0FF,
                                         0x0F00000000F000F0, 0x0F000F0000000F00,
                                         0x000000F000000000, 0x000F00000000F000,
                                         0x00000000000F0000, 0x0000000000000000,
                                         0x0000000000000000]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Certain, a[i]) == certain_a_answers[i]
                        @test Var.create_nibble_mask(Certain, b[i]) == certain_b_answers[i]
                        @test Var.create_nibble_mask(Certain, a[i], b[i]) == (certain_a_answers[i] & certain_b_answers[i])
                    end

                    match_answers = [0x000000000000000F, 0xF00000000000000F,
                                     0x00F0000000000000, 0x00000F0000000000,
                                     0x000000000F000000, 0x0000F000000000F0,
                                     0x000F0000000F0000, 0x0F0000F00000F000,
                                     0xFFFFFFFFF0F00F00]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Match, a[i], b[i]) == match_answers[i]
                    end

                    mismatch_answers = [~i for i in match_answers]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Mismatch, a[i], b[i]) == mismatch_answers[i]
                    end

                    conserved_answers = [0x0000000000000000, 0xF00000000000000F,
                                         0x0000000000000000, 0x00000F0000000000,
                                         0x0000000000000000, 0x0000000000000000,
                                         0x00000000000F0000, 0x0000000000000000,
                                         0x0000000000000000]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Conserved, a[i], b[i]) == conserved_answers[i]
                    end

                    mutated_answers = [0x0000000000000000, 0x00000000F000F0F0,
                                       0x0000000000F000F0, 0x0F00000000000000,
                                       0x0000000000000000, 0x0000000000000000,
                                       0x0000000000000000, 0x0000000000000000,
                                       0x0000000000000000]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Mutated, a[i], b[i]) == mutated_answers[i]
                    end

                    transition_answers = [0x0000000000000000, 0x000000000000F000,
                                          0x0000000000F00000, 0x0000000000000000,
                                          0x0000000000000000, 0x0000000000000000,
                                          0x0000000000000000, 0x0000000000000000,
                                          0x0000000000000000]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Transition, a[i], b[i]) == transition_answers[i]
                    end

                    transversion_answers = [0x0000000000000000, 0x00000000F00000F0,
                                            0x00000000000000F0, 0x0F00000000000000,
                                            0x0000000000000000, 0x0000000000000000,
                                            0x0000000000000000, 0x0000000000000000,
                                            0x0000000000000000]

                    for i in 1:endof(a)
                        @test Var.create_nibble_mask(Transversion, a[i], b[i]) == transversion_answers[i]
                    end
                end

                @testset "Nibble counting" begin

                    @testset "Indels" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Indel, s)
                            off = 16 - n
                            @test (Var.count_sites4(Indel, s.data[1]) - off) == expected
                        end
                    end

                    @testset "Ambiguities" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected_a = count_sites_naieve(Ambiguous, s)
                            expected_b = count_sites_naieve(Ambiguous, s, s2)
                            @test Var.count_sites4(Ambiguous, s.data[1]) == expected_a
                            @test Var.count_sites4(Ambiguous, s.data[1], s2.data[1]) == expected_b
                        end
                    end

                    @testset "Certain" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Certain, s, s2)
                            @test Var.count_sites4(Certain, s.data[1], s2.data[1]) == expected
                        end
                    end

                    @testset "Match" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Match, s, s2)
                            out = 16 - n
                            @test (Var.count_sites4(Match, s.data[1], s2.data[1]) - out) == expected
                        end
                    end

                    @testset "Mismatch" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Mismatch, s, s2)
                            @test Var.count_sites4(Mismatch, s.data[1], s2.data[1]) == expected
                        end
                    end

                    @testset "Conserved" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Conserved, s, s2)
                            @test Var.count_sites4(Conserved, s.data[1], s2.data[1]) == expected
                        end
                    end

                    @testset "Mutated" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Mutated, s, s2)
                            @test Var.count_sites4(Mutated, s.data[1], s2.data[1]) == expected
                        end
                    end

                    @testset "Transition" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Transition, s, s2)
                            @test Var.count_sites4(Transition, s.data[1], s2.data[1]) == expected
                        end
                    end

                    @testset "Transversion" begin
                        for i in 1:500
                            n = rand(1:16)
                            s = generate_testcase(DNAAlphabet{4}, n)
                            s2 = generate_testcase(DNAAlphabet{4}, n)
                            expected = count_sites_naieve(Transversion, s, s2)
                            @test Var.count_sites4(Transversion, s.data[1], s2.data[1]) == expected
                        end
                    end
                end
            end

            @testset "ShiftedIntsItr" begin
                @testset "4 bit biological sequences" begin

                    # Replace test cases with automatic generation of any case.


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
                    # 1200 Mutations, 400 conserved sites, 400 non-certain sites.
                    seq1 = dna"ATTG-ACCTGGNTTTCCGAA"^100
                    seq2 = dna"A-ACAGAGTATACRGTCGTC"^100

                    @testset "Verifying naieve methods" begin
                        @test count_sites_naieve(Indel, seq1, seq2) == 200
                        @test count_sites_naieve(Ambiguous, seq1, seq2) == 200
                        @test count_sites_naieve(Certain, seq1, seq2) == 1600
                        @test count_sites_naieve(Conserved, seq1, seq2) == 400
                        @test count_sites_naieve(Mutated, seq1, seq2) == 1200
                        @test count_sites_naieve(Transition, seq1, seq2) == 400
                        @test count_sites_naieve(Transversion, seq1, seq2) == 800
                    end

                    @testset "Vefifying count_sites4 methods for multiple integers" begin

                        @testset "No shifting" begin
                            # Set cases we've reasoned about to test certain
                            # behaviours and properties.
                            @test Var.count_sites4(Ambiguous, seq1.data, seq2.data) == 200
                            @test Var.count_sites4(Certain, seq1.data, seq2.data) == 1600
                            @test Var.count_sites4(Conserved, seq1.data, seq2.data) == 400
                            @test Var.count_sites4(Mutated, seq1.data, seq2.data) == 1200
                            @test Var.count_sites4(Transition, seq1.data, seq2.data) == 400
                            @test Var.count_sites4(Transversion, seq1.data, seq2.data) == 800
                            # Randomly generated testcases.
                            for i in 1:100
                                n = rand(1:5000)
                                s = generate_testcase(DNAAlphabet{4}, n)
                                s2 = generate_testcase(DNAAlphabet{4}, n)

                                expected_gap = count_sites_naieve(Indel, s, s2)
                                expected_amb = count_sites_naieve(Ambiguous, s, s2)
                                expected_ctn = count_sites_naieve(Certain, s, s2)
                                expected_cns = count_sites_naieve(Conserved, s, s2)
                                expected_mm = count_sites_naieve(Mismatch, s, s2)
                                expected_mut = count_sites_naieve(Mutated, s, s2)
                                expected_trs = count_sites_naieve(Transition, s, s2)
                                expected_trv = count_sites_naieve(Transversion, s, s2)

                                @test Var.count_sites4(Ambiguous, s.data, s2.data) == expected_amb
                                @test Var.count_sites4(Certain, s.data, s2.data) == expected_ctn
                                @test Var.count_sites4(Conserved, s.data, s2.data) == expected_cns
                                @test Var.count_sites4(Mismatch, s.data, s2.data) == expected_mm
                                @test Var.count_sites4(Mutated, s.data, s2.data) == expected_mut
                                @test Var.count_sites4(Transition, s.data, s2.data) == expected_trs
                                @test Var.count_sites4(Transversion, s.data, s2.data) == expected_trv
                            end
                        end



                    end



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
