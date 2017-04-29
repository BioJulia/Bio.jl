module TestAlign

using Base.Test

using Bio
using Bio.Seq
using Bio.Align
using Bio.Intervals
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

        submat = SubstitutionMatrix(DNA, rand(Float64, 15, 15))
        @test isa(submat, SubstitutionMatrix{DNA,Float64})

        submat = DichotomousSubstitutionMatrix(5, -4)
        @test isa(submat, DichotomousSubstitutionMatrix{Int})
        submat = convert(SubstitutionMatrix{DNA,Int}, submat)
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
        submat = SubstitutionMatrix(DNA, rand(Float64, 15, 15))
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
        submat = SubstitutionMatrix(DNA, rand(Int, 15, 15))
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
        samdir = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "SAM")

        @testset "MetaInfo" begin
            metainfo = SAM.MetaInfo()
            @test !isfilled(metainfo)
            @test contains(repr(metainfo), "not filled")

            metainfo = SAM.MetaInfo("CO", "some comment (parens)")
            @test isfilled(metainfo)
            @test string(metainfo) == "@CO\tsome comment (parens)"
            @test contains(repr(metainfo), "CO")
            @test SAM.tag(metainfo) == "CO"
            @test SAM.value(metainfo) == "some comment (parens)"
            @test_throws ArgumentError keys(metainfo)
            @test_throws ArgumentError values(metainfo)

            metainfo = SAM.MetaInfo("HD", ["VN" => "1.0", "SO" => "coordinate"])
            @test isfilled(metainfo)
            @test string(metainfo) == "@HD\tVN:1.0\tSO:coordinate"
            @test contains(repr(metainfo), "HD")
            @test SAM.tag(metainfo) == "HD"
            @test SAM.value(metainfo) == "VN:1.0\tSO:coordinate"
            @test keys(metainfo) == ["VN", "SO"]
            @test values(metainfo) == ["1.0", "coordinate"]
            @test SAM.keyvalues(metainfo) == ["VN" => "1.0", "SO" => "coordinate"]
            @test haskey(metainfo, "VN")
            @test haskey(metainfo, "SO")
            @test !haskey(metainfo, "GO")
            @test metainfo["VN"] == "1.0"
            @test metainfo["SO"] == "coordinate"
            @test_throws KeyError metainfo["GO"]
        end

        @testset "Header" begin
            header = SAM.Header()
            @test isempty(header)
            push!(header, SAM.MetaInfo("@HD\tVN:1.0\tSO:coordinate"))
            @test !isempty(header)
            @test length(header) == 1
            push!(header, SAM.MetaInfo("@CO\tsome comment"))
            @test length(header) == 2
            @test isa(collect(header), Vector{SAM.MetaInfo})
        end

        @testset "Record" begin
            record = SAM.Record()
            @test !isfilled(record)
            @test !SAM.ismapped(record)
            @test repr(record) == "Bio.Align.SAM.Record: <not filled>"
            @test_throws ArgumentError SAM.flag(record)

            record = SAM.Record("r001\t99\tchr1\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*")
            @test isfilled(record)
            @test ismatch(r"^Bio.Align.SAM.Record:\n", repr(record))
            @test SAM.ismapped(record)
            @test SAM.hastempname(record)
            @test SAM.tempname(record) == "r001"
            @test SAM.hasflag(record)
            @test SAM.flag(record) === UInt16(99)
            @test SAM.hasrefname(record)
            @test SAM.refname(record) == "chr1"
            @test SAM.hasposition(record)
            @test SAM.position(record) === 7
            @test SAM.hasmappingquality(record)
            @test SAM.mappingquality(record) === UInt8(30)
            @test SAM.hascigar(record)
            @test SAM.cigar(record) == "8M2I4M1D3M"
            @test SAM.hasnextrefname(record)
            @test SAM.nextrefname(record) == "="
            @test SAM.hasnextposition(record)
            @test SAM.nextposition(record) === 37
            @test SAM.hastemplength(record)
            @test SAM.templength(record) === 39
            @test SAM.hassequence(record)
            @test SAM.sequence(record) == dna"TTAGATAAAGGATACTG"
            @test !SAM.hasquality(record)
            @test_throws MissingFieldException SAM.quality(record)
        end

        @testset "Reader" begin
            reader = open(SAM.Reader, joinpath(samdir, "ce#1.sam"))
            @test isa(reader, SAM.Reader)
            @test eltype(reader) === SAM.Record

            # header
            h = header(reader)
            @test string.(find(header(reader), "SQ")) == ["@SQ\tSN:CHROMOSOME_I\tLN:1009800"]

            # first record
            record = SAM.Record()
            read!(reader, record)
            @test SAM.ismapped(record)
            @test SAM.refname(record) == "CHROMOSOME_I"
            @test SAM.position(record) == leftposition(record) == 2
            @test SAM.rightposition(record) == rightposition(record) == 102
            @test SAM.tempname(record) == seqname(record) == "SRR065390.14978392"
            @test SAM.sequence(record) == sequence(record) == dna"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"
            @test SAM.sequence(String, record)          ==    "CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"
            @test SAM.seqlength(record) == 100
            @test SAM.quality(record)         == (b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" .- 33)
            @test SAM.quality(String, record) ==   "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            @test SAM.flag(record) == 16
            @test SAM.cigar(record) == "27M1D73M"
            @test SAM.alignment(record) == Alignment([
                AlignmentAnchor(  0,   1, OP_START),
                AlignmentAnchor( 27,  28, OP_MATCH),
                AlignmentAnchor( 27,  29, OP_DELETE),
                AlignmentAnchor(100, 102, OP_MATCH)])
            @test record["XG"] == 1
            @test record["XM"] == 5
            @test record["XN"] == 0
            @test record["XO"] == 1
            @test record["AS"] == -18
            @test record["XS"] == -18
            @test record["YT"] == "UU"
            @test eof(reader)
            close(reader)

            # iterator
            @test length(collect(open(SAM.Reader, joinpath(samdir, "ce#1.sam")))) == 1
            @test length(collect(open(SAM.Reader, joinpath(samdir, "ce#2.sam")))) == 2

            # IOStream
            @test length(collect(SAM.Reader(open(joinpath(samdir, "ce#1.sam"))))) == 1
            @test length(collect(SAM.Reader(open(joinpath(samdir, "ce#2.sam"))))) == 2
        end

        @testset "Round trip" begin
            function compare_records(xs, ys)
                if length(xs) != length(ys)
                    return false
                end
                for (x, y) in zip(xs, ys)
                    if x.data[x.filled] != y.data[y.filled]
                        return false
                    end
                end
                return true
            end
            for specimen in YAML.load_file(joinpath(samdir, "index.yml"))
                filepath = joinpath(samdir, specimen["filename"])
                mktemp() do path, io
                    # copy
                    reader = open(SAM.Reader, filepath)
                    writer = SAM.Writer(io, header(reader))
                    records = SAM.Record[]
                    for record in reader
                        push!(records, record)
                        write(writer, record)
                    end
                    close(reader)
                    close(writer)
                    @test compare_records(open(collect, SAM.Reader, path), records)
                end
            end
        end
    end

    @testset "BAM" begin
        bamdir = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "BAM")

        @testset "AuxData" begin
            auxdata = BAM.AuxData(UInt8[])
            @test isempty(auxdata)

            buf = IOBuffer()
            write(buf, "NM", UInt8('s'), Int16(1))
            auxdata = BAM.AuxData(takebuf_array(buf))
            @test length(auxdata) == 1
            @test auxdata["NM"] === Int16(1)
            @test collect(auxdata) == ["NM" => Int16(1)]

            buf = IOBuffer()
            write(buf, "AS", UInt8('c'), Int8(-18))
            write(buf, "NM", UInt8('s'), Int16(1))
            write(buf, "XA", UInt8('f'), Float32(3.14))
            write(buf, "XB", UInt8('Z'), "some text\0")
            write(buf, "XC", UInt8('B'), UInt8('i'), Int32(3), Int32[10, -5, 8])
            auxdata = BAM.AuxData(takebuf_array(buf))
            @test length(auxdata) == 5
            @test auxdata["AS"] === Int8(-18)
            @test auxdata["NM"] === Int16(1)
            @test auxdata["XA"] === Float32(3.14)
            @test auxdata["XB"] == "some text"
            @test auxdata["XC"] == Int32[10, -5, 8]
            @test convert(Dict{String,Any}, auxdata) == Dict(
                "AS" => Int8(-18),
                "NM" => Int16(1),
                "XA" => Float32(3.14),
                "XB" => "some text",
                "XC" => Int32[10, -5, 8])
        end

        @testset "Record" begin
            record = BAM.Record()
            @test !isfilled(record)
            @test repr(record) == "Bio.Align.BAM.Record: <not filled>"
            @test_throws ArgumentError BAM.flag(record)
        end

        @testset "Reader" begin
            reader = open(BAM.Reader, joinpath(bamdir, "ce#1.bam"))
            @test isa(reader, BAM.Reader)
            @test eltype(reader) === BAM.Record
            @test startswith(repr(reader), "Bio.Align.BAM.Reader{IOStream}:")

            # header
            h = header(reader)
            @test isa(h, SAM.Header)

            # first record
            record = BAM.Record()
            read!(reader, record)
            @test BAM.ismapped(record)
            @test BAM.refname(record) == "CHROMOSOME_I"
            @test BAM.refid(record) === 1
            @test BAM.hasnextrefid(record)
            @test BAM.nextrefid(record) === 0
            @test BAM.hasposition(record) === hasleftposition(record) === true
            @test BAM.position(record) === leftposition(record) === 2
            @test BAM.hasnextposition(record)
            @test BAM.nextposition(record) === 0
            @test rightposition(record) == 102
            @test BAM.hastempname(record) === hasseqname(record) === true
            @test BAM.tempname(record) == seqname(record) == "SRR065390.14978392"
            @test BAM.hassequence(record) === hassequence(record) === true
            @test BAM.sequence(record) == sequence(record) == dna"""
            CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCT
            AAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
            """
            @test BAM.seqlength(record) === 100
            @test BAM.hasquality(record)
            @test eltype(BAM.quality(record)) == UInt8
            @test BAM.quality(record) == [Int(x) - 33 for x in "#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"]
            @test BAM.flag(record) === UInt16(16)
            @test BAM.cigar(record) == "27M1D73M"
            @test BAM.alignment(record) == Alignment([
                AlignmentAnchor(  0,   1, OP_START),
                AlignmentAnchor( 27,  28, OP_MATCH),
                AlignmentAnchor( 27,  29, OP_DELETE),
                AlignmentAnchor(100, 102, OP_MATCH)])
            @test record["XG"] == 1
            @test record["XM"] == 5
            @test record["XN"] == 0
            @test record["XO"] == 1
            @test record["AS"] == -18
            @test record["XS"] == -18
            @test record["YT"] == "UU"
            @test keys(record) == ["XG","XM","XN","XO","AS","XS","YT"]
            @test values(record) == [1, 5, 0, 1, -18, -18, "UU"]
            @test eof(reader)
            close(reader)

            # iterator
            @test length(collect(open(BAM.Reader, joinpath(bamdir, "ce#1.bam")))) == 1
            @test length(collect(open(BAM.Reader, joinpath(bamdir, "ce#2.bam")))) == 2

            # IOStream
            @test length(collect(BAM.Reader(open(joinpath(bamdir, "ce#1.bam"))))) == 1
            @test length(collect(BAM.Reader(open(joinpath(bamdir, "ce#2.bam"))))) == 2
        end

        function compare_records(xs, ys)
            if length(xs) != length(ys)
                return false
            end
            for (x, y) in zip(xs, ys)
                if !(
                    x.block_size == y.block_size &&
                    x.refid      == y.refid &&
                    x.pos        == y.pos &&
                    x.bin_mq_nl  == y.bin_mq_nl &&
                    x.flag_nc    == y.flag_nc &&
                    x.l_seq      == y.l_seq &&
                    x.next_refid == y.next_refid &&
                    x.next_pos   == y.next_pos &&
                    x.tlen       == y.tlen &&
                    x.data[1:BAM.data_size(x)] == y.data[1:BAM.data_size(y)])
                    return false
                end
            end
            return true
        end

        @testset "Round trip" begin
            for specimen in YAML.load_file(joinpath(bamdir, "index.yml"))
                filepath = joinpath(bamdir, specimen["filename"])
                mktemp() do path, _
                    # copy
                    if contains(get(specimen, "tags", ""), "bai")
                        reader = open(BAM.Reader, filepath, index=filepath * ".bai")
                    else
                        reader = open(BAM.Reader, filepath)
                    end
                    writer = BAM.Writer(
                        BGZFStream(path, "w"),
                        BAM.header(reader, fillSQ=isempty(find(header(reader), "SQ"))))
                    records = BAM.Record[]
                    for record in reader
                        push!(records, record)
                        write(writer, record)
                    end
                    close(reader)
                    close(writer)
                    @test compare_records(open(collect, BAM.Reader, path), records)
                end
            end
        end

        @testset "Random access" begin
            filepath = joinpath(bamdir, "GSE25840_GSM424320_GM06985_gencode_spliced.head.bam")
            reader = open(BAM.Reader, filepath, index=filepath * ".bai")

            @test isa(eachoverlap(reader, "chr1", 1:100), BAM.OverlapIterator)
            @test isa(eachoverlap(reader, Interval("chr1", 1, 100)), BAM.OverlapIterator)

            # expected values are counted using samtools
            for (refname, interval, expected) in [
                    ("chr1", 1_000:10000,      21),
                    ("chr1", 8_000:10000,      20),
                    ("chr1", 766_000:800_000, 142),
                    ("chr1", 786_000:800_000, 1),
                    ("chr1", 796_000:800_000, 0)]
                intsect = eachoverlap(reader, refname, interval)
                @test eltype(intsect) == BAM.Record
                @test count(_ -> true, intsect) == expected
                # check that the intersection iterator is stateless
                @test count(_ -> true, intsect) == expected
            end

            # randomized tests
            for n in 1:50
                refindex = 1
                refname = "chr1"
                range = randrange(1:1_000_000)
                seekstart(reader)
                # linear scan
                expected = collect(filter(reader) do record
                    Bio.Align.BAM.isoverlapping(record, refindex, range)
                end)
                # indexed scan
                actual = collect(eachoverlap(reader, refname, range))
                @test compare_records(actual, expected)
            end
            close(reader)

            filepath = joinpath(bamdir, "R_12h_D06.uniq.q40.bam")
            reader = open(BAM.Reader, filepath, index=filepath * ".bai")
            @test isempty(collect(eachoverlap(reader, "chr19", 5823708:5846478)))
            close(reader)
        end
    end
end

end # TestAlign
