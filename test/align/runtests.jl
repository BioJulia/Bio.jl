module TestAlign

using FactCheck
using Bio
using Bio.Seq
using Bio.Align
using TestFunctions


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


facts("Alignments") do
    context("Operations") do
        context("Constructors and Conversions") do
            ops = Set(Align.char_to_op)
            for op in ops
                if op != Align.OP_INVALID
                    @fact Operation(Char(op)) --> op
                end
            end
            @fact_throws Char(Align.OP_INVALID)
            @fact_throws Operation('m')
            @fact_throws Operation('7')
            @fact_throws Operation('A')
            @fact_throws Operation('\n')
        end
    end

    context("AlignmentAnchor") do
        anchor = AlignmentAnchor(1, 2, OP_MATCH)
        @fact string(anchor) --> "AlignmentAnchor(1, 2, 'M')"
    end

    context("Alignment") do
        # alignments with nonsense operations
        @fact_throws Alignment(AlignmentAnchor[
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
            @fact_throws Alignment(anchors)
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
            @fact_throws Alignment(anchors)
        end

        # cigar string round-trip
        for _ in 1:100
            path = random_alignment(rand(1000:10000), rand(1000:10000))
            anchors = anchors_from_path(path)
            aln = Alignment(anchors)
            cig = cigar(aln)
            @fact Alignment(cig, aln.anchors[1].seqpos + 1,
                            aln.anchors[1].refpos + 1) --> aln
        end
    end

    context("AlignedSequence") do
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
        @fact Bio.Align.first(alnseq) -->  5
        @fact Bio.Align.last(alnseq)  --> 27
        # OP_MATCH
        for (seqpos, refpos) in [(1, 5), (2, 6), (4, 8), (13, 18), (19, 27)]
            @fact seq2ref(seqpos, alnseq) --> (refpos, OP_MATCH)
            @fact ref2seq(refpos, alnseq) --> (seqpos, OP_MATCH)
        end
        # OP_INSERT
        @fact seq2ref(10, alnseq) --> (17, OP_INSERT)
        @fact seq2ref(11, alnseq) --> (17, OP_INSERT)
        # OP_DELETE
        @fact ref2seq( 9, alnseq) --> ( 4, OP_DELETE)
        @fact ref2seq(10, alnseq) --> ( 4, OP_DELETE)
        @fact ref2seq(23, alnseq) --> (15, OP_DELETE)
    end
end


# generate test cases from two aligned sequences
function alnscore{S,T}(::Type{S}, affinegap::AffineGapScoreModel{T}, alnstr::ASCIIString, clip::Bool)
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

function alnscore{T}(affinegap::AffineGapScoreModel{T}, alnstr::ASCIIString; clip=true)
    return alnscore(ASCIIString, affinegap, alnstr, clip)
end

function alndistance{S,T}(::Type{S}, cost::CostModel{T}, alnstr::ASCIIString)
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

function alndistance{T}(cost::CostModel{T}, alnstr::ASCIIString)
    return alndistance(ASCIIString, cost, alnstr)
end

function alignedpair(alnres)
    a, b = alignment(alnres)
    anchors = a.aln.anchors
    buf = IOBuffer()
    Bio.Align.show_seq(buf, a, anchors)
    println(buf)
    Bio.Align.show_ref(buf, b, anchors)
    return bytestring(buf)
end

facts("PairwiseAlignment") do
    context("SubstitutionMatrix") do
        # defined
        @fact BLOSUM62[AA_A,AA_R] --> -1
        @fact BLOSUM62[AA_R,AA_A] --> -1
        @fact BLOSUM62[AA_R,AA_R] -->  5
        @fact typeof(BLOSUM62[AA_A,AA_R]) --> Int
        # undefined
        @fact BLOSUM62[AA_O,AA_R] -->  0
        @fact BLOSUM62[AA_R,AA_O] -->  0
        @fact Bio.Align.is_defined_symbol(BLOSUM62, AA_R) --> true
        @fact Bio.Align.is_defined_symbol(BLOSUM62, AA_O) --> false

        # no error
        print(IOBuffer(), BLOSUM62)

        mat = Int[
            +2 -3 -3 -3  0;
            -3 +2 -3 -3  0;
            -3 -3 +2 -3  0;
            -3 -3 -3 +2  0;
             0  0  0  0  0;
        ]
        defined = BitVector([1, 1, 1, 1, 0])
        match = 0
        mismatch = -1
        submat = SubstitutionMatrix{Int}(mat, defined, match, mismatch)
        @fact submat[DNA_A,DNA_A] --> +2
        @fact submat[DNA_N,DNA_N] -->  0
        @fact submat[DNA_N,DNA_A] --> -1
        @fact submat[DNA_A,DNA_N] --> -1
    end

    context("AffineGapScoreModel") do
        # predefined substitution matrix
        for affinegap in [AffineGapScoreModel(BLOSUM62, -10, -1),
                          AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1),
                          AffineGapScoreModel(BLOSUM62, gap_open_penalty=10, gap_extend_penalty=1)]
            @fact affinegap.gap_open --> -10
            @fact affinegap.gap_extend --> -1
            @fact typeof(affinegap) --> AffineGapScoreModel{Int}
        end

        # matrix
        submat = Float64[
             1 -1 -1 -1;
            -1  1 -1 -1;
            -1 -1  1 -1;
            -1 -1 -1  1;
        ]
        for affinegap in [AffineGapScoreModel(submat, -3, -1),
                          AffineGapScoreModel(submat, gap_open=-3, gap_extend=-1),
                          AffineGapScoreModel(submat, gap_open_penalty=3, gap_extend_penalty=1)]
            @fact affinegap.gap_open --> -3
            @fact affinegap.gap_extend --> -1
            @fact typeof(affinegap) --> AffineGapScoreModel{Float64}
        end

        affinegap = AffineGapScoreModel(match=3, mismatch=-3, gap_open=-5, gap_extend=-2)
        @fact affinegap.gap_open --> -5
        @fact affinegap.gap_extend --> -2
        @fact typeof(affinegap) --> AffineGapScoreModel{Int}
    end

    context("CostModel") do
        submat = [
            0 3 3 3;
            3 0 3 3;
            3 3 0 3;
            3 3 3 0;
        ]
        for cost in [CostModel(submat, 5, 6),
                     CostModel(submat, insertion=5, deletion=6)]
            @fact cost.insertion --> 5
            @fact cost.deletion --> 6
            @fact typeof(cost) --> CostModel{Int}
        end

        cost = CostModel(match=0, mismatch=3, insertion=5, deletion=6)
        @fact cost.insertion --> 5
        @fact cost.deletion --> 6
        @fact typeof(cost) --> CostModel{Int}
    end

    context("Alignment") do
        anchors = [
            AlignmentAnchor(0, 0, OP_START),
            AlignmentAnchor(3, 3, OP_SEQ_MATCH)
        ]
        seq = AlignedSequence("ACG", anchors)
        ref = "ACG"
        aln = PairwiseAlignment(seq, ref)
        a, b = aln
        @fact a === seq --> true
        @fact b === ref --> true
        result = PairwiseAlignmentResult(3, true, seq, ref)
        @fact isa(result, PairwiseAlignmentResult) --> true
        @fact isa(alignment(result), PairwiseAlignment) --> true
        @fact score(result) --> 3
        @fact hasalignment(result) --> true
    end

    context("count_<ops>") do
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
        @fact count_matches(aln) --> 7
        @fact count_mismatches(aln) --> 4
        @fact count_insertions(aln) --> 2
        @fact count_deletions(aln) --> 3
        @fact count_aligned(aln) --> 16
    end

    context("GlobalAlignment") do
        affinegap = AffineGapScoreModel(
            match=0,
            mismatch=-6,
            gap_open=-5,
            gap_extend=-3
        )

        function testaln(alnstr)
            a, b, s, alnpair = alnscore(affinegap, alnstr)
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact score(aln) --> s
            @fact alignedpair(aln) --> alnpair
            aln = pairalign(GlobalAlignment(), a, b, affinegap, score_only=true)
            @fact score(aln) --> s
        end

        context("empty sequences") do
            aln = pairalign(GlobalAlignment(), "", "", affinegap)
            @fact score(aln) --> 0
        end

        context("complete match") do
            testaln("""
            ACGT
            ACGT
            """)
        end

        context("mismatch") do
            testaln("""
            ACGT
            AGGT
            """)

            testaln("""
            ACGT
            AGGA
            """)
        end

        context("insertion") do
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

        context("deletion") do
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

        context("banded") do
            a, b, s, alnpair = alnscore(affinegap, """
            ACGT
            ACGT
            """)
            aln = pairalign(GlobalAlignment(), a, b, affinegap, banded=true)
            @fact score(aln) --> s
            @fact alignedpair(aln) --> alnpair

            a, b, s, alnpair = alnscore(affinegap, """
            ACGT
            AGGT
            """)
            aln = pairalign(GlobalAlignment(), a, b, affinegap, banded=true)
            @fact score(aln) --> s
            @fact alignedpair(aln) --> alnpair

            a, b, s, alnpair = alnscore(affinegap, """
            ACG--T
            ACGAAT
            """)
            aln = pairalign(GlobalAlignment(), a, b, affinegap, banded=true, lower_offset=0, upper_offset=0)
            @fact score(aln) --> s
            @fact alignedpair(aln) --> alnpair
        end
    end

    context("SemiGlobalAlignment") do
        affinegap = AffineGapScoreModel(
            match=0,
            mismatch=-6,
            gap_open=-5,
            gap_extend=-3
        )

        function testaln(alnstr)
            a, b, s, alnpair = alnscore(affinegap, alnstr, clip=false)
            aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
            @fact score(aln) --> s
            @fact alignedpair(aln) --> alnpair
            aln = pairalign(SemiGlobalAlignment(), a, b, affinegap, score_only=true)
            @fact score(aln) --> s
        end

        context("complete match") do
            testaln("""
            ACGT
            ACGT
            """)
        end

        context("partial match") do
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

    context("OverlapAlignment") do
        affinegap = AffineGapScoreModel(
            match=3,
            mismatch=-6,
            gap_open=-5,
            gap_extend=-3
        )

        function testaln(alnstr)
            a, b, s, alnpair = alnscore(affinegap, alnstr, clip=false)
            aln = pairalign(OverlapAlignment(), a, b, affinegap)
            @fact score(aln) --> s
            @fact alignedpair(aln) --> alnpair
            aln = pairalign(OverlapAlignment(), a, b, affinegap, score_only=true)
            @fact score(aln) --> s
        end

        context("complete match") do
            testaln("""
            ACGT
            ACGT
            """)
        end

        context("partial match") do
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

    context("LocalAlignment") do
        context("zero matching score") do
            affinegap = AffineGapScoreModel(
                match=0,
                mismatch=-6,
                gap_open=-5,
                gap_extend=-3
            )

            function testaln(alnstr)
                a, b, s, alnpair = alnscore(affinegap, alnstr)
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact score(aln) --> s
                @fact alignedpair(aln) --> alnpair
                aln = pairalign(LocalAlignment(), a, b, affinegap, score_only=true)
                @fact score(aln) --> s
            end

            context("empty sequences") do
                aln = pairalign(LocalAlignment(), "", "", affinegap)
                @fact score(aln) --> 0
            end

            context("complete match") do
                testaln("""
                ACGT
                ACGT
                """)
            end

            context("partial match") do
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

            context("no match") do
                a = "AA"
                b = "TTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact score(aln) --> 0
            end
        end

        context("positive matching score") do
            affinegap = AffineGapScoreModel(
                match=5,
                mismatch=-6,
                gap_open=-5,
                gap_extend=-3
            )

            function testaln(alnstr)
                a, b, s, alnpair = alnscore(affinegap, alnstr)
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact score(aln) --> s
                @fact alignedpair(aln) --> alnpair
                aln = pairalign(LocalAlignment(), a, b, affinegap, score_only=true)
                @fact score(aln) --> s
            end

            context("complete match") do
                testaln("""
                ACGT
                ACGT
                ^^^^
                """)
            end

            context("partial match") do
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

            context("no match") do
                a = "AA"
                b = "TTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact score(aln) --> 0
            end
        end
    end

    context("EditDistance") do
        mismatch = 1
        submat = DichotomousSubstitutionMatrix(0, mismatch)
        insertion = 1
        deletion = 2
        cost = CostModel(submat, insertion, deletion)

        function testaln(alnstr)
            a, b, dist = alndistance(cost, alnstr)
            aln = pairalign(EditDistance(), a, b, cost)
            @fact distance(aln) --> dist
            @fact alignedpair(aln) --> chomp(alnstr)
            aln = pairalign(EditDistance(), a, b, cost, distance_only=true)
            @fact distance(aln) --> dist
        end

        context("empty sequences") do
            aln = pairalign(EditDistance(), "", "", cost)
            @fact distance(aln) --> 0
        end

        context("complete match") do
            testaln("""
            ACGT
            ACGT
            """)
        end

        context("mismatch") do
            testaln("""
            AGGT
            ACGT
            """)

            testaln("""
            AGGT
            ACGT
            """)
        end

        context("insertion") do
            testaln("""
            ACGTT
            ACG-T
            """)

            testaln("""
            ACGTT
            -CG-T
            """)
        end

        context("deletion") do
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

    context("LevenshteinDistance") do
        context("empty sequences") do
            aln = pairalign(LevenshteinDistance(), "", "")
            @fact distance(aln) --> 0
        end

        context("complete match") do
            a = "ACGT"
            b = "ACGT"
            aln = pairalign(LevenshteinDistance(), a, b)
            @fact distance(aln) --> 0
        end
    end

    context("HammingDistance") do
        function testaln(alnstr)
            a, b = split(chomp(alnstr), '\n')
            dist = sum([x != y for (x, y) in zip(a, b)])
            aln = pairalign(HammingDistance(), a, b)
            @fact distance(aln) --> dist
            @fact alignedpair(aln) --> chomp(alnstr)
            aln = pairalign(HammingDistance(), a, b, distance_only=true)
            @fact distance(aln) --> dist
        end

        context("empty sequences") do
            aln = pairalign(HammingDistance(), "", "")
            @fact distance(aln) --> 0
        end

        context("complete match") do
            testaln("""
            ACGT
            ACGT
            """)
        end

        context("mismatch") do
            testaln("""
            ACGT
            AGGT
            """)

            testaln("""
            ACGT
            AGGA
            """)
        end

        context("indel") do
            @fact_throws pairalign(HammingDistance(), "ACGT", "ACG")
            @fact_throws pairalign(HammingDistance(), "ACG", "ACGT")
        end
    end
end

end # TestAlign
