# Pairwise-Alignment Result
# -------------------------

typealias PairedSequences{S1,S2} Pair{AlignedSequence{S1},S2}

type PairwiseAlignment{T,S1,S2}
    score::T
    seqpair::Nullable{PairedSequences{S1,S2}}
end

function Base.call{S1,S2}(::Type{PairwiseAlignment{S1,S2}}, score)
    return PairwiseAlignment(score, Nullable{PairedSequences{S1,S2}}())
end

function PairwiseAlignment(score, seq, ref)
    return PairwiseAlignment(score, Nullable(seq => ref))
end

function Base.show{T,S1,S2}(io::IO, aln::PairwiseAlignment{T,S1,S2})
    print(io, "PairwiseAlignment{", T, ",", S1, ",", S2, "}:", '\n')
    print(io, "  score: ", aln.score)
    if !isnull(aln.seqpair)
        pair = get(aln.seqpair)
        seq = pair.first
        ref = pair.second
        print(io, '\n')
        anchors = seq.aln.anchors
        # show the aligned sequence
        print(io, "  seq: ")
        for i in 2:length(anchors)
            if ismatchop(anchors[i].op) || isinsertop(anchors[i].op)
                for j in anchors[i-1].seqpos+1:anchors[i].seqpos
                    print(io, seq.seq[j])
                end
            elseif isdeleteop(anchors[i].op)
                for _ in anchors[i-1].refpos+1:anchors[i].refpos
                    write(io, '-')
                end
            end
        end
        println(io)
        # show the matching string
        print(io, "       ")
        for i in 2:length(anchors)
            op = anchors[i].op
            if ismatchop(op)
                for _ in anchors[i-1].seqpos+1:anchors[i].seqpos
                    print(io, op == OP_SEQ_MATCH ? '|' : ' ')
                end
            elseif isinsertop(op)
                for _ in anchors[i-1].seqpos+1:anchors[i].seqpos
                    print(io, ' ')
                end
            elseif isdeleteop(op)
                for _ in anchors[i-1].refpos+1:anchors[i].refpos
                    print(io, ' ')
                end
            end
        end
        println(io)
        # show the reference sequence
        print(io, "  ref: ")
        for i in 2:length(anchors)
            if ismatchop(anchors[i].op) || isdeleteop(anchors[i].op)
                for j in anchors[i-1].refpos+1:anchors[i].refpos
                    print(io, ref[j])
                end
            elseif isinsertop(anchors[i].op)
                for _ in anchors[i-1].seqpos+1:anchors[i].seqpos
                    write(io, '-')
                end
            end
        end
    end
end
