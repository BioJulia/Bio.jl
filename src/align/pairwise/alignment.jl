# Pairwise Alignment
# ------------------

"""
Pairwise alignment
"""
type PairwiseAlignment{S1,S2}
    a::AlignedSequence{S1}
    b::S2
end

Base.start(aln::PairwiseAlignment) = 1
Base.done(aln::PairwiseAlignment, i) = i > 2
Base.next(aln::PairwiseAlignment, i) = (i == 1 ? aln.a : aln.b), i + 1

"""
The number of matching positions.
"""
function nmatching(aln::PairwiseAlignment)
    seq = aln.a
    anchors = seq.aln.anchors
    n = 0
    for i in 2:endof(anchors)
        op = anchors[i].op
        if op == OP_SEQ_MATCH
            n += anchors[i].seqpos - anchors[i-1].seqpos
        end
    end
    return n
end

"""
The number of mismatching positions.
"""
function nmismatching(aln::PairwiseAlignment)
    seq = aln.a
    anchors = seq.aln.anchors
    n = 0
    for i in 2:endof(anchors)
        op = anchors[i].op
        if ismatchop(op) && op != OP_SEQ_MATCH
            n += anchors[i].seqpos - anchors[i-1].seqpos
        end
    end
    return n
end

"""
The number of aligned positions.
"""
function naligned(aln::PairwiseAlignment)
    seq = aln.a
    anchors = seq.aln.anchors
    n = 0
    for i in 2:endof(anchors)
        op = anchors[i].op
        if ismatchop(op) || isinsertop(op)
            n += anchors[i].seqpos - anchors[i-1].seqpos
        elseif isdeleteop(op)
            n += anchors[i].refpos - anchors[i-1].refpos
        end
    end
    return n
end

function Base.show{S1,S2}(io::IO, aln::PairwiseAlignment{S1,S2})
    println(io, "PairwiseAlignment{", S1, ",", S2, "}:")
    show_pairwise_alignment(io, aln)
end

function show_pairwise_alignment(io::IO, aln::PairwiseAlignment)
    seq, ref = aln
    anchors = seq.aln.anchors
    println(io)
    print(io, "  seq: "); show_seq(io, seq, anchors); println(io)
    print(io, "       "); show_match(io, anchors); println(io)
    print(io, "  ref: "); show_ref(io, ref, anchors)
end

function show_seq(io, seq, anchors)
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
end

function show_match(io, anchors)
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
end

function show_ref(io, ref, anchors)
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
