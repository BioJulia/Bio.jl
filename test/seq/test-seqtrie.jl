@testset "SeqTrie" begin
    seqs = Dict(
        DNAAlphabet{4} => [dna"ACGT", dna"ACGTACGT", dna"TGCA"],
        RNAAlphabet{4} => [rna"ACGU", rna"ACGUACGU", rna"UGCA"],
        AminoAcidAlphabet => [aa"PWFS", aa"PWFSPWFS", aa"WCATDWG"],
    )

    @testset "Constructors" begin
        # Alphabet aliases
        @test typeof(DNATrie()) == SeqTrie{DNAAlphabet{4}, Any}
        @test typeof(RNATrie()) == SeqTrie{RNAAlphabet{4}, Any}
        @test typeof(AminoAcidTrie()) == SeqTrie{AminoAcidAlphabet, Any}
    end

    @testset "Alias types" begin
        d = DNATrie()
        r = RNATrie()
        a = AminoAcidTrie()

        for (i, seqtype) in enumerate((DNASequence, RNASequence, AminoAcidSequence))
            seq = seqtype("ACGA") # can be DNA, RNA or AA
            for (j, trie) in enumerate((d, r, a))
                if i == j
                    trie[seq] = 1
                    @test trie[seq] == 1
                else
                    @test_throws Exception trie[seq] = 1
                end
            end
            dict = Dict{seqtype, Int}(seq => 1)
            trie = SeqTrie(dict)
            trie[seq] = 1
            @test trie[seq] == 1
        end
    end

    @testset "Basic Operations" for AB in (DNAAlphabet{4}, RNAAlphabet{4}, AminoAcidAlphabet)
        t = SeqTrie{AB, Int}()

        for (i, seq) in enumerate(seqs[AB])
            t[seq] = i
        end

        for (i, seq) in enumerate(seqs[AB])
            @test haskey(t, seq)
            @test get(t, seq, nothing) == i
        end

        @test sort(keys(t)) == seqs[AB]
        @test sort(keys_with_prefix(t, seqs[AB][1])) == seqs[AB][1:2]
    end

    @testset "Prefix iterator" begin
        t = DNATrie()

        for (i, seq) in enumerate(seqs[DNAAlphabet{4}])
            t[seq] = i
        end
        t0 = t
        t1 = t0.children[DNA_A]
        t2 = t1.children[DNA_C]
        t3 = t2.children[DNA_G]
        # empty => root
        @test collect(path(t, dna"")) == [t0]
        # not present => root
        @test collect(path(t, dna"G")) == [t0]
        # prefix => correct nodes
        @test collect(path(t, dna"ACG")) == [t0, t1, t2, t3]
        # prefix present, key absent => correct nodes of prefix
        @test collect(path(t, dna"ACGG")) == [t0, t1, t2, t3]
    end
end
