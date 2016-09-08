module TestAlign

using Base.Test

using Bio
using Bio.Seq
using Bio.Align
using BGZFStreams
using TestFunctions
using YAML


# Generate a random valid alignment of a sequence of length n against a sequence
# of length m. If `glob` is true, generate a global alignment, if false, a local
# alignment.
function random_alignment(m, n, glob=true)
    match_ops = [OP_MATCH, OP_SEQ_MATCH, OP_SEQ_MISMATCH]
    insert_ops = [OP_INSERT, OP_SOFT_CLIP, OP_HARD_CLIP]
    delete_ops = [OP_DELETE, OP_SKIP]
    ops = vcat(match_ops, insert_ops, delete_ops)

    # This is just a random walk on a m-by-n matrix, where steps are either
    # (+1,0), (0,+1), (+1,+1). To make somewhat more realistic alignments, it's
    # biased towards going in the same direction. Local alignments have a random
    # start and end time, global alignments always start at (0,0) and end at
    # (m,n).

    # probability of choosing the same direction as the last step
    straight_pr = 0.9

    op = OP_MATCH
    if glob
        i = 0
        j = 0
        i_end = m
        j_end = n
    else
        i = rand(1:m-1)
        j = rand(1:n-1)
        i_end = rand(i+1:m)
        j_end = rand(j+1:n)
    end

    path = AlignmentAnchor[AlignmentAnchor(i, j, OP_START)]
    while (glob && i < i_end && j < j_end) || (!glob && (i < i_end || j < j_end))
        straight = rand() < straight_pr

        if i == i_end
            if !straight
                op = rand(delete_ops)
            end
            j += 1
        elseif j == j_end
            if !straight
                op = rand(inset_ops)
            end
            i += 1
        else
            if !straight
                op = rand(ops)
            end

            if isdeleteop(op)
                j += 1
            elseif isinsertop(op)
                i += 1
            else
                i += 1
                j += 1
            end
        end
        push!(path, AlignmentAnchor(i, j, op))
    end

    return path
end


# Make an Alignment from a path returned by random_alignment. Converting from
# path to Alignment is just done by removing redundant nodes from the path.
function anchors_from_path(path)
    anchors = AlignmentAnchor[]
    for k in 1:length(path)
        if k == length(path) || path[k].op != path[k+1].op
            push!(anchors, path[k])
        end
    end
    return anchors
end

# Generate a random range within `range`.
function randrange(range)
    x = rand(range)
    y = rand(range)
    if x < y
        return x:y
    else
        return y:x
    end
end

@testset "Alignments" begin
    @testset "Operations" begin
        @testset "Constructors and Conversions" begin
            ops = Set(Align.char_to_op)
            for op in ops
                if op != Align.OP_INVALID
                    @test Operation(Char(op)) == op
                end
            end
            @test_throws Exception Char(Align.OP_INVALID)
            @test_throws Exception Operation('m')
            @test_throws Exception Operation('7')
            @test_throws Exception Operation('A')
            @test_throws Exception Operation('\n')
        end
    end

    @testset "AlignmentAnchor" begin
        anchor = AlignmentAnchor(1, 2, OP_MATCH)
        @test string(anchor) == "AlignmentAnchor(1, 2, 'M')"
    end

    @testset "Alignment" begin
        # alignments with nonsense operations
        @test_throws Exception Alignment(AlignmentAnchor[
            Operation(0, 0, OP_START),
            Operation(100, 100, convert(Operation, 0xfa))])

        # test bad alignment anchors by swapping nodes in paths
        for _ in 1:100
            path = random_alignment(rand(1000:10000), rand(1000:10000))
            anchors = anchors_from_path(path)
            n = length(anchors)
            n < 3 && continue
            i = rand(2:n-1)
            j = rand(i+1:n)
            anchors[i], anchors[j] = anchors[j], anchors[i]
            @test_throws Exception Alignment(anchors)
        end

        # test bad alignment anchors by swapping operations
        for _ in 1:100
            path = random_alignment(rand(1000:10000), rand(1000:10000))
            anchors = anchors_from_path(path)
            n = length(anchors)
            n < 3 && continue
            i = rand(2:n-1)
            j = rand(i+1:n)
            u = anchors[i]
            v = anchors[j]
            if (ismatchop(u.op) && ismatchop(v.op)) ||
               (isinsertop(u.op) && isinsertop(v.op)) ||
               (isdeleteop(u.op) && isdeleteop(v.op))
                continue
            end
            anchors[i] = AlignmentAnchor(u.seqpos, u.refpos, v.op)
            anchors[j] = AlignmentAnchor(v.seqpos, v.refpos, u.op)
            @test_throws Exception Alignment(anchors)
        end

        # cigar string round-trip
        for _ in 1:100
            path = random_alignment(rand(1000:10000), rand(1000:10000))
            anchors = anchors_from_path(path)
            aln = Alignment(anchors)
            cig = cigar(aln)
            @test Alignment(cig, aln.anchors[1].seqpos + 1,
                            aln.anchors[1].refpos + 1) == aln
        end
    end

    @testset "AlignedSequence" begin
        #               0   4        9  12 15     19
        #               |   |        |  |  |      |
        #     query:     TGGC----ATCATTTAACG---CAAG
        # reference: AGGGTGGCATTTATCAG---ACGTTTCGAGAC
        #               |   |   |    |     |  |   |
        #               4   8   12   17    20 23  27
        anchors = [
            AlignmentAnchor( 0,  4, OP_START),
            AlignmentAnchor( 4,  8, OP_MATCH),
            AlignmentAnchor( 4, 12, OP_DELETE),
            AlignmentAnchor( 9, 17, OP_MATCH),
            AlignmentAnchor(12, 17, OP_INSERT),
            AlignmentAnchor(15, 20, OP_MATCH),
            AlignmentAnchor(15, 23, OP_DELETE),
            AlignmentAnchor(19, 27, OP_MATCH)
        ]
        query = "TGGCATCATTTAACGCAAG"
        alnseq = AlignedSequence(query, anchors)
        @test Bio.Align.first(alnseq) ==  5
        @test Bio.Align.last(alnseq)  == 27
        # OP_MATCH
        for (seqpos, refpos) in [(1, 5), (2, 6), (4, 8), (13, 18), (19, 27)]
            @test seq2ref(alnseq, seqpos) == (refpos, OP_MATCH)
            @test ref2seq(alnseq, refpos) == (seqpos, OP_MATCH)
        end
        # OP_INSERT
        @test seq2ref(alnseq, 10) == (17, OP_INSERT)
        @test seq2ref(alnseq, 11) == (17, OP_INSERT)
        # OP_DELETE
        @test ref2seq(alnseq,  9) == ( 4, OP_DELETE)
        @test ref2seq(alnseq, 10) == ( 4, OP_DELETE)
        @test ref2seq(alnseq, 23) == (15, OP_DELETE)

        seq = dna"ACGG--TGAAAGGT"
        ref = dna"-CGGGGA----TTT"
        alnseq = AlignedSequence(seq, ref)
        @test Bio.Align.first(alnseq) == 1
        @test Bio.Align.last(alnseq)  == 9
        @test alnseq.aln.anchors == [
             AlignmentAnchor( 0, 0, '0')
             AlignmentAnchor( 1, 0, 'I')
             AlignmentAnchor( 4, 3, '=')
             AlignmentAnchor( 4, 5, 'D')
             AlignmentAnchor( 5, 6, 'X')
             AlignmentAnchor( 9, 6, 'I')
             AlignmentAnchor(11, 8, 'X')
             AlignmentAnchor(12, 9, '=')
        ]
    end
end


# generate test cases from two aligned sequences
function alnscore{S,T}(::Type{S}, affinegap::AffineGapScoreModel{T}, alnstr::AbstractString, clip::Bool)
    gap_open = affinegap.gap_open
    gap_extend = affinegap.gap_extend
    lines = split(chomp(alnstr), '\n')
    a, b = lines[1:2]
    m = length(a)
    @assert m == length(b)

    if length(lines) == 2
        start = 1
        while start ≤ m && a[start] == ' ' || b[start] == ' '
            start += 1
        end
        stop = start
        while stop + 1 ≤ m && !(a[stop+1] == ' ' || b[stop+1] == ' ')
            stop += 1
        end
    elseif length(lines) == 3
        start = search(lines[3], '^')
        stop = rsearch(lines[3], '^')
    else
        error("invalid alignment string")
    end

    score = T(0)
    gap_extending_a = false
    gap_extending_b = false
    for i in start:stop
        if a[i] == '-'
            score += gap_extending_a ? gap_extend : (gap_open + gap_extend)
            gap_extending_a = true
        elseif b[i] == '-'
            score += gap_extending_b ? gap_extend : (gap_open + gap_extend)
            gap_extending_b = true
        else
            score += affinegap.submat[a[i],b[i]]
            gap_extending_a = false
            gap_extending_b = false
        end
    end
    sa = S(replace(a, r"\s|-", ""))
    sb = S(replace(b, r"\s|-", ""))
    return sa, sb, score, clip ? string(a[start:stop], '\n', b[start:stop]) : string(a, '\n', b)
end

function alnscore{T}(affinegap::AffineGapScoreModel{T}, alnstr::AbstractString; clip=true)
    return alnscore(AbstractString, affinegap, alnstr, clip)
end

function alndistance{S,T}(::Type{S}, cost::CostModel{T}, alnstr::AbstractString)
    lines = split(chomp(alnstr), '\n')
    @assert length(lines) == 2
    a, b = lines
    m = length(a)
    @assert length(b) == m
    dist = T(0)
    for i in 1:m
        if a[i] == '-'
            dist += cost.deletion
        elseif b[i] == '-'
            dist += cost.insertion
        else
            dist += cost.submat[a[i],b[i]]
        end
    end
    return S(replace(a, r"\s|-", "")), S(replace(b, r"\s|-", "")), dist
end

function alndistance{T}(cost::CostModel{T}, alnstr::AbstractString)
    return alndistance(AbstractString, cost, alnstr)
end

function alignedpair(alnres)
    aln = alignment(alnres)
    a = aln.a
    b = aln.b
    anchors = a.aln.anchors
    buf = IOBuffer()
    print_seq(buf, a, anchors)
    println(buf)
    print_ref(buf, b, anchors)
    return String(buf)
end

function print_seq(io, seq, anchors)
    for i in 2:length(anchors)
        if ismatchop(anchors[i].op) || isinsertop(anchors[i].op)
            for j in anchors[i-1].seqpos+1:anchors[i].seqpos
                print(io, seq.seq[j])
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                print(io, '-')
            end
        end
    end
end

function print_ref(io, ref, anchors)
    for i in 2:length(anchors)
        if ismatchop(anchors[i].op) || isdeleteop(anchors[i].op)
            for j in anchors[i-1].refpos+1:anchors[i].refpos
                print(io, ref[j])
            end
        elseif isinsertop(anchors[i].op)
            for _ in anchors[i-1].seqpos+1:anchors[i].seqpos
                print(io, '-')
            end
        end
    end
end

@testset "PairwiseAlignment" begin
    @testset "SubstitutionMatrix" begin
        # DNA
        @test EDNAFULL[DNA_A,DNA_A] ===  5
        @test EDNAFULL[DNA_G,DNA_G] ===  5
        @test EDNAFULL[DNA_A,DNA_G] === -4
        @test EDNAFULL[DNA_G,DNA_A] === -4
        @test EDNAFULL[DNA_M,DNA_T] === -4
        @test EDNAFULL[DNA_M,DNA_C] ===  1

        # amino acid
        @test BLOSUM62[AA_A,AA_R] === -1
        @test BLOSUM62[AA_R,AA_A] === -1
        @test BLOSUM62[AA_R,AA_R] ===  5
        @test BLOSUM62[AA_O,AA_R] ===  0  # default
        @test BLOSUM62[AA_R,AA_O] ===  0  # default

        # update
        myblosum = copy(BLOSUM62)
        @test myblosum[AA_A,AA_R] === -1
        myblosum[AA_A,AA_R] = 10
        @test myblosum[AA_A,AA_R] === 10

        @test BLOSUM62[AA_O,AA_R] ===  0  # default
        myblosum[AA_O,AA_R] = -3
        @test myblosum[AA_O,AA_R] === -3

        submat = SubstitutionMatrix(DNANucleotide, rand(Float64, 15, 15))
        @test isa(submat, SubstitutionMatrix{DNANucleotide,Float64})

        submat = DichotomousSubstitutionMatrix(5, -4)
        @test isa(submat, DichotomousSubstitutionMatrix{Int})
        submat = convert(SubstitutionMatrix{DNANucleotide,Int}, submat)
        @test submat[DNA_A,DNA_A] ===  5
        @test submat[DNA_C,DNA_C] ===  5
        @test submat[DNA_A,DNA_C] === -4
        @test submat[DNA_C,DNA_A] === -4

        try
            print(IOBuffer(), EDNAFULL)
            print(IOBuffer(), BLOSUM62)
            # no error
            @test true
        catch
            @test false
        end
    end

    @testset "AffineGapScoreModel" begin
        # predefined substitution matrix
        for affinegap in [AffineGapScoreModel(BLOSUM62, -10, -1),
                          AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1),
                          AffineGapScoreModel(BLOSUM62, gap_open_penalty=10, gap_extend_penalty=1)]
            @test affinegap.gap_open == -10
            @test affinegap.gap_extend == -1
            @test typeof(affinegap) == AffineGapScoreModel{Int}
        end

        # matrix
        submat = SubstitutionMatrix(DNANucleotide, rand(Float64, 15, 15))
        for affinegap in [AffineGapScoreModel(submat, -3, -1),
                          AffineGapScoreModel(submat, gap_open=-3, gap_extend=-1),
                          AffineGapScoreModel(submat, gap_open_penalty=3, gap_extend_penalty=1)]
            @test affinegap.gap_open == -3
            @test affinegap.gap_extend == -1
            @test typeof(affinegap) == AffineGapScoreModel{Float64}
        end

        affinegap = AffineGapScoreModel(match=3, mismatch=-3, gap_open=-5, gap_extend=-2)
        @test affinegap.gap_open == -5
        @test affinegap.gap_extend == -2
        @test typeof(affinegap) == AffineGapScoreModel{Int}
    end

    @testset "CostModel" begin
        submat = SubstitutionMatrix(DNANucleotide, rand(Int, 15, 15))
        for cost in [CostModel(submat, 5, 6),
                     CostModel(submat, insertion=5, deletion=6)]
            @test cost.insertion == 5
            @test cost.deletion == 6
            @test typeof(cost) == CostModel{Int}
        end

        cost = CostModel(match=0, mismatch=3, insertion=5, deletion=6)
        @test cost.insertion == 5
        @test cost.deletion == 6
        @test typeof(cost) == CostModel{Int}
    end

    @testset "Alignment" begin
        anchors = [
            AlignmentAnchor(0, 0, OP_START),
            AlignmentAnchor(3, 3, OP_SEQ_MATCH)
        ]
        seq = AlignedSequence("ACG", anchors)
        ref = "ACG"
        aln = PairwiseAlignment(seq, ref)
        @test collect(aln) == [('A', 'A'), ('C', 'C'), ('G', 'G')]
        result = PairwiseAlignmentResult(3, true, seq, ref)
        @test isa(result, PairwiseAlignmentResult) == true
        @test isa(alignment(result), PairwiseAlignment) == true
        @test score(result) == 3
        @test hasalignment(result) == true
    end

    @testset "count_<ops>" begin
        # anchors are derived from an alignment:
        #   seq: ACG---TGCAGAATTT
        #        |     || || ||  
        #   ref: AAAATTTGAAGTAT--
        a = dna"ACGTGCAGAATTT"
        b = dna"AAAATTTGAAGTAT"
        anchors = [
            AlignmentAnchor( 0,  0, '0'),
            AlignmentAnchor( 1,  1, '='),
            AlignmentAnchor( 3,  3, 'X'),
            AlignmentAnchor( 3,  6, 'D'),
            AlignmentAnchor( 5,  8, '='),
            AlignmentAnchor( 6,  9, 'X'),
            AlignmentAnchor( 8, 11, '='),
            AlignmentAnchor( 9, 12, 'X'),
            AlignmentAnchor(11, 14, '='),
            AlignmentAnchor(13, 14, 'I')
        ]
        aln = PairwiseAlignment(AlignedSequence(a, anchors), b)
        @test count_matches(aln) == 7
        @test count_mismatches(aln) == 4
        @test count_insertions(aln) == 2
        @test count_deletions(aln) == 3
        @test count_aligned(aln) == 16
    end

    @testset "Interfaces" begin
        seq = dna"ACGTATAGT"
        ref = dna"ATCGTATTGGT"
        # seq:  1 A-CGTATAG-T  9
        #         | ||||| | |
        # ref:  1 ATCGTATTGGT 11
        model = AffineGapScoreModel(EDNAFULL, gap_open=-4, gap_extend=-1)
        result = pairalign(GlobalAlignment(), seq, ref, model)
        @test isa(result, PairwiseAlignmentResult)
        aln = alignment(result)
        @test isa(aln, PairwiseAlignment)
        @test seq2ref(aln, 1) == (1, OP_SEQ_MATCH)
        @test seq2ref(aln, 2) == (3, OP_SEQ_MATCH)
        @test seq2ref(aln, 3) == (4, OP_SEQ_MATCH)
        @test ref2seq(aln, 1) == (1, OP_SEQ_MATCH)
        @test ref2seq(aln, 2) == (1, OP_DELETE)
        @test ref2seq(aln, 3) == (2, OP_SEQ_MATCH)
    end

    @testset "GlobalAlignment" begin
        affinegap = AffineGapScoreModel(
            match=0,
            mismatch=-6,
            gap_open=-5,
            gap_extend=-3
        )

        function testaln(alnstr)
            a, b, s, alnpair = alnscore(affinegap, alnstr)
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @test score(aln) == s
            @test alignedpair(aln) == alnpair
            aln = pairalign(GlobalAlignment(), a, b, affinegap, score_only=true)
            @test score(aln) == s
        end

        @testset "empty sequences" begin
            aln = pairalign(GlobalAlignment(), "", "", affinegap)
            @test score(aln) == 0
        end

        @testset "complete match" begin
            testaln("""
            ACGT
            ACGT
            """)
        end

        @testset "mismatch" begin
            testaln("""
            ACGT
            AGGT
            """)

            testaln("""
            ACGT
            AGGA
            """)
        end

        @testset "insertion" begin
            testaln("""
            ACGTT
            ACGT-
            """)

            testaln("""
            ACGTTT
            ACGT--
            """)

            testaln("""
            ACCGT
            AC-GT
            """)

            testaln("""
            ACCCGT
            AC--GT
            """)

            testaln("""
            AACGT
            A-CGT
            """)

            testaln("""
            AAACGT
            A--CGT
            """)
        end

        @testset "deletion" begin
            testaln("""
            ACGT-
            ACGTT
            """)

            testaln("""
            ACGT-
            ACGTT
            """)

            testaln("""
            ACGT--
            ACGTTT
            """)

            testaln("""
            AC-GT
            ACCGT
            """)

            testaln("""
            AC--GT
            ACCCGT
            """)

            testaln("""
            A-CGT
            AACGT
            """)

            testaln("""
            A--CGT
            AAACGT
            """)
        end

        @testset "banded" begin
            a, b, s, alnpair = alnscore(affinegap, """
            ACGT
            ACGT
            """)
            aln = pairalign(GlobalAlignment(), a, b, affinegap, banded=true)
            @test score(aln) == s
            @test alignedpair(aln) == alnpair

            a, b, s, alnpair = alnscore(affinegap, """
            ACGT
            AGGT
            """)
            aln = pairalign(GlobalAlignment(), a, b, affinegap, banded=true)
            @test score(aln) == s
            @test alignedpair(aln) == alnpair

            a, b, s, alnpair = alnscore(affinegap, """
            ACG--T
            ACGAAT
            """)
            aln = pairalign(GlobalAlignment(), a, b, affinegap, banded=true, lower_offset=0, upper_offset=0)
            @test score(aln) == s
            @test alignedpair(aln) == alnpair
        end
    end

    @testset "SemiGlobalAlignment" begin
        affinegap = AffineGapScoreModel(
            match=0,
            mismatch=-6,
            gap_open=-5,
            gap_extend=-3
        )

        function testaln(alnstr)
            a, b, s, alnpair = alnscore(affinegap, alnstr, clip=false)
            aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
            @test score(aln) == s
            @test alignedpair(aln) == alnpair
            aln = pairalign(SemiGlobalAlignment(), a, b, affinegap, score_only=true)
            @test score(aln) == s
        end

        @testset "complete match" begin
            testaln("""
            ACGT
            ACGT
            """)
        end

        @testset "partial match" begin
            testaln("""
            --ACTT---
            TTACGTAGT
              ^^^^
            """)

            testaln("""
            --AC-TTG-
            TTACGTTGT
              ^^^^^^
            """)

            testaln("""
            --ACTAGT---
            TTAC--GTTGT
              ^^^^^^
            """)
        end
    end

    @testset "OverlapAlignment" begin
        affinegap = AffineGapScoreModel(
            match=3,
            mismatch=-6,
            gap_open=-5,
            gap_extend=-3
        )

        function testaln(alnstr)
            a, b, s, alnpair = alnscore(affinegap, alnstr, clip=false)
            aln = pairalign(OverlapAlignment(), a, b, affinegap)
            @test score(aln) == s
            @test alignedpair(aln) == alnpair
            aln = pairalign(OverlapAlignment(), a, b, affinegap, score_only=true)
            @test score(aln) == s
        end

        @testset "complete match" begin
            testaln("""
            ACGT
            ACGT
            """)
        end

        @testset "partial match" begin
            testaln("""
            ---ACGGTGATTAT
            GATACGGTGA----
               ^^^^^^^
            """)

            testaln("""
            ---AACGT-GATTAT
            GATAACGGAGA----
               ^^^^^^^^
            """)

            testaln("""
            GATACGGTGA----
            ---ACGGTGATTAT
               ^^^^^^^
            """)

            testaln("""
            GATAACGGAGA----
            ---AACGT-GATTAT
               ^^^^^^^^
            """)
        end
    end

    @testset "LocalAlignment" begin
        @testset "zero matching score" begin
            affinegap = AffineGapScoreModel(
                match=0,
                mismatch=-6,
                gap_open=-5,
                gap_extend=-3
            )

            function testaln(alnstr)
                a, b, s, alnpair = alnscore(affinegap, alnstr)
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @test score(aln) == s
                @test alignedpair(aln) == alnpair
                aln = pairalign(LocalAlignment(), a, b, affinegap, score_only=true)
                @test score(aln) == s
            end

            @testset "empty sequences" begin
                aln = pairalign(LocalAlignment(), "", "", affinegap)
                @test score(aln) == 0
            end

            @testset "complete match" begin
                testaln("""
                ACGT
                ACGT
                """)
            end

            @testset "partial match" begin
                testaln("""
                ACGT
                AGGT
                  ^^
                """)

                testaln("""
                   ACGT
                AACGTTT
                      ^
                """)
            end

            @testset "no match" begin
                a = "AA"
                b = "TTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @test score(aln) == 0
            end
        end

        @testset "positive matching score" begin
            affinegap = AffineGapScoreModel(
                match=5,
                mismatch=-6,
                gap_open=-5,
                gap_extend=-3
            )

            function testaln(alnstr)
                a, b, s, alnpair = alnscore(affinegap, alnstr)
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @test score(aln) == s
                @test alignedpair(aln) == alnpair
                aln = pairalign(LocalAlignment(), a, b, affinegap, score_only=true)
                @test score(aln) == s
            end

            @testset "complete match" begin
                testaln("""
                ACGT
                ACGT
                ^^^^
                """)
            end

            @testset "partial match" begin
                testaln("""
                ACGT
                AGGT
                  ^^
                """)

                testaln("""
                 ACGT  
                AACGTTT
                 ^^^^
                """)

                testaln("""
                  AC-GT  
                AAACTGTTT
                """)
            end

            @testset "no match" begin
                a = "AA"
                b = "TTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @test score(aln) == 0
            end
        end
    end

    @testset "EditDistance" begin
        mismatch = 1
        submat = DichotomousSubstitutionMatrix(0, mismatch)
        insertion = 1
        deletion = 2
        cost = CostModel(submat, insertion, deletion)

        function testaln(alnstr)
            a, b, dist = alndistance(cost, alnstr)
            aln = pairalign(EditDistance(), a, b, cost)
            @test distance(aln) == dist
            @test alignedpair(aln) == chomp(alnstr)
            aln = pairalign(EditDistance(), a, b, cost, distance_only=true)
            @test distance(aln) == dist
        end

        @testset "empty sequences" begin
            aln = pairalign(EditDistance(), "", "", cost)
            @test distance(aln) == 0
        end

        @testset "complete match" begin
            testaln("""
            ACGT
            ACGT
            """)
        end

        @testset "mismatch" begin
            testaln("""
            AGGT
            ACGT
            """)

            testaln("""
            AGGT
            ACGT
            """)
        end

        @testset "insertion" begin
            testaln("""
            ACGTT
            ACG-T
            """)

            testaln("""
            ACGTT
            -CG-T
            """)
        end

        @testset "deletion" begin
            testaln("""
            AC-T
            ACGT
            """)

            testaln("""
            -C-T
            ACGT
            """)
        end
    end

    @testset "LevenshteinDistance" begin
        @testset "empty sequences" begin
            aln = pairalign(LevenshteinDistance(), "", "")
            @test distance(aln) == 0
        end

        @testset "complete match" begin
            a = "ACGT"
            b = "ACGT"
            aln = pairalign(LevenshteinDistance(), a, b)
            @test distance(aln) == 0
        end
    end

    @testset "HammingDistance" begin
        function testaln(alnstr)
            a, b = split(chomp(alnstr), '\n')
            dist = sum([x != y for (x, y) in zip(a, b)])
            aln = pairalign(HammingDistance(), a, b)
            @test distance(aln) == dist
            @test alignedpair(aln) == chomp(alnstr)
            aln = pairalign(HammingDistance(), a, b, distance_only=true)
            @test distance(aln) == dist
        end

        @testset "empty sequences" begin
            aln = pairalign(HammingDistance(), "", "")
            @test distance(aln) == 0
        end

        @testset "complete match" begin
            testaln("""
            ACGT
            ACGT
            """)
        end

        @testset "mismatch" begin
            testaln("""
            ACGT
            AGGT
            """)

            testaln("""
            ACGT
            AGGA
            """)
        end

        @testset "indel" begin
            @test_throws Exception pairalign(HammingDistance(), "ACGT", "ACG")
            @test_throws Exception pairalign(HammingDistance(), "ACG", "ACGT")
        end
    end
end

@testset "High-throughput Sequencing" begin
    get_bio_fmt_specimens()

    @testset "AuxDataDict" begin
        dict = AuxDataDict()
        @test length(dict) == 0
        @test isempty(dict)
        @test_throws KeyError dict["NM"]

        dict = AuxDataDict(
            "X1" => UInt8(1),
            "X2" => UInt16(2),
            "X3" => UInt32(3),
            "X4" => Int8(4),
            "X5" => Int16(5),
            "X6" => Int32(6),
            "X7" => Float32(7),
            "X8" => "eight",
            "X9" => Int32[9])
        @test length(dict) == 9
        @test !isempty(dict)
        @test dict["X1"] === UInt8(1)
        @test dict["X2"] === UInt16(2)
        @test dict["X3"] === UInt32(3)
        @test dict["X4"] === Int8(4)
        @test dict["X5"] === Int16(5)
        @test dict["X6"] === Int32(6)
        @test dict["X7"] === Float32(7)
        @test dict["X8"] == "eight"
        @test typeof(dict["X8"]) == String
        @test dict["X9"] == Int32[9]
        @test typeof(dict["X9"]) == Vector{Int32}

        dict = AuxDataDict("NM" => 0x01, "XY" => Int32(100), "XZ" => [0x11, 0x23])
        @test length(dict) == 3
        @test dict["NM"] === 0x01
        @test dict["XY"] === Int32(100)
        @test dict["XZ"] == [0x11, 0x23]
        @test eltype(dict["XZ"]) == UInt8

        dict = AuxDataDict("NM" => 0x01, "MD" => "8T1T39")
        @test length(dict) == 2
        @test dict["NM"] === 0x01
        @test dict["MD"] == "8T1T39"
        dict["NM"] = 0x00
        @test dict["NM"] === 0x00
        dict["MD"] = "50"
        @test dict["MD"] == "50"
        @test collect(dict) == ["NM" => 0x00, "MD" => "50"]
        dict["XY"] = "foobar"
        @test dict["XY"] == "foobar"
        @test collect(dict) == ["NM" => 0x00, "MD" => "50", "XY" => "foobar"]
        delete!(dict, "NM")
        @test length(dict) == 2
        @test collect(dict) == ["MD" => "50", "XY" => "foobar"]
    end

    @testset "SAM" begin
        samdir = Pkg.dir("Bio", "test", "BioFmtSpecimens", "SAM")

        @testset "SAMHeader" begin
            h = SAMHeader()
            @test isa(h, Associative)
            @test isempty(h)
            h["HD"] = Dict("VN" => "100.100", "SO" => "unknown")
            @test length(h) == 1
            @test h["HD"]["VN"] == "100.100"
            h["CO"] = ["comment1", "comment2"]
            @test length(h) == 2
            @test h["CO"] == ["comment1", "comment2"]
            delete!(h, "CO")
            @test length(h) == 1
        end

        @testset "Record" begin
            rec = SAMRecord()
            @test !ismapped(rec)

            # default values
            @test seqname(rec) == "*"
            @test flag(rec) == SAM_FLAG_UNMAP
            @test refname(rec) == "*"
            @test leftposition(rec) == 0
            @test rightposition(rec) == -1
            @test nextrefname(rec) == "*"
            @test nextleftposition(rec) == 0
            @test templatelength(rec) == 0
            @test cigar(rec) == "*"
            @test alignment(rec) == Alignment(AlignmentAnchor[])
            @test sequence(rec) == "*"
            @test qualities(rec) == "*"
            @test Align.alignment_length(rec) === 0

            # set & delete tags
            rec = SAMRecord()
            @test !haskey(rec, "MN")
            rec["MN"] = 0x01
            @test rec["MN"] === 0x01
            @test haskey(rec, "MN")
            @test !haskey(rec, "XY")
            rec["XY"] = "foobar"
            @test rec["XY"] == "foobar"
            @test haskey(rec, "XY")
            delete!(rec, "MN")
            @test !haskey(rec, "MN")
        end

        @testset "Reader" begin
            reader = open(SAMReader, joinpath(samdir, "ce#1.sam"))
            @test isa(reader, SAMReader)
            @test eltype(reader) === SAMRecord
            @test startswith(repr(reader), "Bio.Align.SAMReader:")

            # header
            h = header(reader)
            @test h["SQ"] == [Dict("SN" => "CHROMOSOME_I", "LN" => "1009800")]

            # first record
            rec = SAMRecord()
            read!(reader, rec)
            @test ismapped(rec)
            @test refname(rec) == "CHROMOSOME_I"
            @test leftposition(rec) == 2
            @test rightposition(rec) == 102
            @test seqname(rec) == "SRR065390.14978392"
            @test sequence(rec)  == "CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"
            @test qualities(rec) == "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            @test flag(rec) == 16
            @test cigar(rec) == "27M1D73M"
            @test alignment(rec) == Alignment([
                AlignmentAnchor(  0,   1, OP_START),
                AlignmentAnchor( 27,  28, OP_MATCH),
                AlignmentAnchor( 27,  29, OP_DELETE),
                AlignmentAnchor(100, 102, OP_MATCH)])
            @test rec["XG"] == 1
            @test rec["XM"] == 5
            @test rec["XN"] == 0
            @test rec["XO"] == 1
            @test rec["AS"] == -18
            @test rec["XS"] == -18
            @test rec["YT"] == "UU"
            @test eof(reader)
            close(reader)

            # iterator
            @test length(collect(open(SAMReader, joinpath(samdir, "ce#1.sam")))) == 1
            @test length(collect(open(SAMReader, joinpath(samdir, "ce#2.sam")))) == 2

            # IOStream
            @test length(collect(SAMReader(open(joinpath(samdir, "ce#1.sam"))))) == 1
            @test length(collect(SAMReader(open(joinpath(samdir, "ce#2.sam"))))) == 2
        end

        @testset "Round trip" begin
            for specimen in YAML.load_file(joinpath(samdir, "index.yml"))
                filepath = joinpath(samdir, specimen["filename"])
                mktemp() do path, io
                    # copy
                    reader = open(SAMReader, filepath)
                    writer = SAMWriter(io, header(reader))
                    records = SAMRecord[]
                    for rec in reader
                        push!(records, rec)
                        write(writer, rec)
                    end
                    close(reader)
                    close(writer)

                    # read again
                    reader = open(SAMReader, path)
                    @test collect(reader) == records
                    close(reader)
                end
            end
        end
    end

    @testset "BAM" begin
        bamdir = Pkg.dir("Bio", "test", "BioFmtSpecimens", "BAM")

        @testset "Record" begin
            rec = BAMRecord()
            @test !ismapped(rec)

            # default values
            @test refname(rec) == "*"
            @test refindex(rec) == 0
            @test leftposition(rec) == 0
            @test rightposition(rec) == -1
            @test mappingquality(rec) == 0
            @test flag(rec) == SAM_FLAG_UNMAP
            @test nextrefname(rec) == "*"
            @test nextrefindex(rec) == 0
            @test nextleftposition(rec) == 0
            @test templatelength(rec) == 0
            @test seqname(rec) == ""
            @test cigar(rec) == ""
            @test sequence(rec) == dna""
            @test alignment(rec) == Alignment(AlignmentAnchor[])
            @test qualities(rec) == UInt8[]
            @test Align.alignment_length(rec) === 0

            # set & delete tags
            rec = BAMRecord()
            @test !haskey(rec, "MN")
            rec["MN"] = 0x01
            @test rec["MN"] === 0x01
            @test haskey(rec, "MN")
            @test !haskey(rec, "XY")
            rec["XY"] = "foobar"
            @test rec["XY"] == "foobar"
            @test haskey(rec, "XY")
            delete!(rec, "MN")
            @test !haskey(rec, "MN")
        end

        @testset "Reader" begin
            reader = open(BAMReader, joinpath(bamdir, "ce#1.bam"))
            @test isa(reader, BAMReader)
            @test eltype(reader) === BAMRecord
            @test startswith(repr(reader), "Bio.Align.BAMReader:")

            # header
            h = header(reader)
            @test h["SQ"] == [Dict("SN" => "CHROMOSOME_I", "LN" => "1009800")]

            # first record
            rec = BAMRecord()
            read!(reader, rec)
            @test ismapped(rec)
            @test refname(rec) == "CHROMOSOME_I"
            @test refindex(rec) == 1
            @test leftposition(rec) == 2
            @test rightposition(rec) == 102
            @test seqname(rec) == "SRR065390.14978392"
            @test sequence(rec) == dna"""
            CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCT
            AAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
            """
            @test eltype(qualities(rec)) == Int8
            @test qualities(rec) == [Int(x) - 33 for x in "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"]
            @test flag(rec) == 16
            @test cigar(rec) == "27M1D73M"
            @test alignment(rec) == Alignment([
                AlignmentAnchor(  0,   1, OP_START),
                AlignmentAnchor( 27,  28, OP_MATCH),
                AlignmentAnchor( 27,  29, OP_DELETE),
                AlignmentAnchor(100, 102, OP_MATCH)])
            @test rec["XG"] == 1
            @test rec["XM"] == 5
            @test rec["XN"] == 0
            @test rec["XO"] == 1
            @test rec["AS"] == -18
            @test rec["XS"] == -18
            @test rec["YT"] == "UU"
            @test eof(reader)
            close(reader)

            # iterator
            @test length(collect(open(BAMReader, joinpath(bamdir, "ce#1.bam")))) == 1
            @test length(collect(open(BAMReader, joinpath(bamdir, "ce#2.bam")))) == 2

            # IOStream
            @test length(collect(BAMReader(open(joinpath(bamdir, "ce#1.bam"))))) == 1
            @test length(collect(BAMReader(open(joinpath(bamdir, "ce#2.bam"))))) == 2
        end

        @testset "Round trip" begin
            for specimen in YAML.load_file(joinpath(bamdir, "index.yml"))
                filepath = joinpath(bamdir, specimen["filename"])
                mktemp() do path, _
                    # copy
                    reader = open(BAMReader, filepath)
                    writer = BAMWriter(
                        BGZFStream(path, "w"),
                        header(reader, true))
                    records = BAMRecord[]
                    for rec in reader
                        push!(records, rec)
                        write(writer, rec)
                    end
                    close(reader)
                    close(writer)

                    # read again
                    reader = open(BAMReader, path)
                    @test collect(reader) == records
                    close(reader)
                end
            end
        end

        @testset "Random access" begin
            filepath = joinpath(bamdir, "GSE25840_GSM424320_GM06985_gencode_spliced.head.bam")
            reader = open(BAMReader, filepath, index=filepath * ".bai")

            # expected values are counted using samtools
            for (refname, interval, expected) in [
                    ("chr1", 1000:10000,      21),
                    ("chr1", 8000:10000,      20),
                    ("chr1", 766_000:800_000, 142),
                    ("chr1", 786_000:800_000, 1),
                    ("chr1", 796_000:800_000, 0)]
                n = 0
                for rec in intersect(reader, refname, interval)
                    n += 1
                end
                @test n == expected
            end

            # randomized tests
            for n in 1:50
                refindex = 1
                refname = "chr1"
                range = randrange(1:1_000_000)
                seekstart(reader)
                # linear scan
                expected = collect(filter(reader) do rec
                    isoverlapping(rec, refindex, range)
                end)
                # indexed scan
                actual = collect(intersect(reader, refname, range))
                @test actual == expected
            end

            close(reader)
        end
    end
end

end # TestAlign
