module TestAlign

using FactCheck
using Bio
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
            AlignmentAnchor(0, 4, OP_START),
            AlignmentAnchor(4, 8, OP_MATCH),
            AlignmentAnchor(4, 12, OP_DELETE),
            AlignmentAnchor(9, 17, OP_MATCH),
            AlignmentAnchor(12, 17, OP_INSERT),
            AlignmentAnchor(15, 20, OP_MATCH),
            AlignmentAnchor(15, 23, OP_DELETE),
            AlignmentAnchor(19, 27, OP_MATCH)
        ]
        query = "TGGCATCATTTAACGCAAG"
        alnseq = AlignedSequence(query, anchors)
        @fact Bio.Align.first(alnseq) --> 5
        @fact Bio.Align.last(alnseq) --> 27
    end
end

# remove gaps
macro g_str(s)
    replace(s, r"-", "")
end

facts("PairwiseAlignment") do
    context("GlobalAlignment") do
        mismatch_score = -6
        gap_open_penalty = 5
        gap_extend_penalty = 3
        submat = DichotomousSubstitutionMatrix(0, mismatch_score)
        affinegap = AffineGapScoreModel(submat, gap_open_penalty, gap_extend_penalty)

        context("empty sequences") do
            a = ""
            b = ""
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> 0
        end

        context("complete match") do
            a = "ACGT"
            b = "ACGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> 0
        end

        context("mismatch") do
            a = "ACGT"
            b = "AGGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> mismatch_score

            a = "ACGT"
            b = "AGGA"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> 2mismatch_score
        end

        context("insertion") do
            a = "ACGTT"
            b = "ACGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + gap_extend_penalty)

            a = "ACGTTT"
            b = "ACGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + 2gap_extend_penalty)

            a = "ACCGT"
            b = "ACGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + gap_extend_penalty)

            a = "ACCCGT"
            b = "ACGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + 2gap_extend_penalty)
        end

        context("deletion") do
            a = "ACGT"
            b = "ACG"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + gap_extend_penalty)

            a = "ACGT"
            b = "AC"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + 2gap_extend_penalty)

            a = "AGT"
            b = "ACGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + gap_extend_penalty)

            a = "AT"
            b = "ACGT"
            aln = pairalign(GlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + 2gap_extend_penalty)
        end
    end

    context("SemiGlobalAlignment") do
        mismatch_score = -6
        gap_open_penalty = 5
        gap_extend_penalty = 3
        submat = DichotomousSubstitutionMatrix(0, mismatch_score)
        affinegap = AffineGapScoreModel(submat, gap_open_penalty, gap_extend_penalty)

        context("compleme match") do
            a = "ACGT"
            b = "ACGT"
            aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
            @fact aln.score --> 0
        end

        context("partial match") do
            a =   "ACTT"
            b = "TTACGTAGT"
            aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
            @fact aln.score --> mismatch_score

            a =  g"AC-TTG"
            b = "TTACGTTGT"
            aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
            @fact aln.score --> -(gap_open_penalty + gap_extend_penalty)
        end
    end

    context("LocalAlignment") do
        context("zero matching score") do
            mismatch_score = -6
            gap_open_penalty = 5
            gap_extend_penalty = 3
            submat = DichotomousSubstitutionMatrix(0, mismatch_score)
            affinegap = AffineGapScoreModel(submat, gap_open_penalty, gap_extend_penalty)

            context("empty sequences") do
                a = ""
                b = ""
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 0
            end

            context("complete match") do
                a = "ACGT"
                b = "ACGT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 0
            end

            context("partial match") do
                a = "ACGT"
                b = "AGGT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 0

                a =  "ACGT"
                b = "AACGTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 0
            end

            context("no match") do
                a = "AA"
                b = "TTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 0
            end
        end

        context("positive matching score") do
            match_score = 5
            mismatch_score = -6
            gap_open_penalty = 5
            gap_extend_penalty = 3
            submat = DichotomousSubstitutionMatrix(match_score, mismatch_score)
            affinegap = AffineGapScoreModel(submat, gap_open_penalty, gap_extend_penalty)

            context("complete match") do
                a = "ACGT"
                b = "ACGT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 4match_score
            end

            context("partial match") do
                a = "ACGT"
                b = "AGGT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 2match_score

                a =  "ACGT"
                b = "AACGTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 4match_score

                a =  g"AC-GT"
                b = "AAACTGTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 4match_score - (gap_open_penalty + gap_extend_penalty)
            end

            context("no match") do
                a = "AA"
                b = "TTTT"
                aln = pairalign(LocalAlignment(), a, b, affinegap)
                @fact aln.score --> 0
            end
        end
    end

    context("EditDistance") do
        mismatch = 1
        submat = DichotomousSubstitutionMatrix(0, mismatch)
        insertion = 1
        deletion = 2
        cost = CostModel(submat, insertion, deletion)

        context("empty sequences") do
            a = ""
            b = ""
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> 0
        end

        context("complete match") do
            a = "ACGT"
            b = "ACGT"
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> 0
        end

        context("mismatch") do
            a = "AGGT"
            b = "ACGT"
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> mismatch

            a = "AGGT"
            b = "ACGA"
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> 2mismatch
        end

        context("insertion") do
            a =  "ACGTT"
            b = g"ACGT-"
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> insertion

            a =  "ACGTT"
            b = g"-CGT-"
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> 2insertion
        end

        context("deletion") do
            a = g"AC-T"
            b =  "ACGT"
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> deletion

            a = g"C-T"
            b = "ACGT"
            aln = pairalign(EditDistance(), a, b, cost)
            @fact aln.score --> 2deletion
        end
    end

    context("LevenshteinDistance") do
        context("empty sequences") do
            a = ""
            b = ""
            aln = pairalign(LevenshteinDistance(), a, b)
            @fact aln.score --> 0
        end

        context("complete match") do
            a = "ACGT"
            b = "ACGT"
            aln = pairalign(LevenshteinDistance(), a, b)
            @fact aln.score --> 0
        end
    end

    context("HammingDistance") do
        context("empty sequences") do
            a = ""
            b = ""
            aln = pairalign(HammingDistance(), a, b)
            @fact aln.score --> 0
        end

        context("complete match") do
            a = "ACGT"
            b = "ACGT"
            aln = pairalign(HammingDistance(), a, b)
            @fact aln.score --> 0
        end

        context("mismatch") do
            a = "ACGT"
            b = "AGGT"
            aln = pairalign(HammingDistance(), a, b)
            @fact aln.score --> 1

            a = "ACGT"
            b = "AGGA"
            aln = pairalign(HammingDistance(), a, b)
            @fact aln.score --> 2
        end

        context("indel") do
            a = "ACGT"
            b = "ACG"
            @fact_throws pairalign(HammingDistance(), a, b)

            a = "ACG"
            b = "ACGT"
            @fact_throws pairalign(HammingDistance(), a, b)
        end
    end
end

end # TestAlign
