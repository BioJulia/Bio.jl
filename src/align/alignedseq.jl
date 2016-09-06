# AlignedSequence
# ===============
#
# An aligned sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable AlignedSequence{S}
    seq::S
    aln::Alignment
end

function AlignedSequence{S}(seq::S, anchors::Vector{AlignmentAnchor},
                            check::Bool=true)
    return AlignedSequence(seq, Alignment(anchors, check))
end

function AlignedSequence(seq::Sequence, ref::Sequence)
    return AlignedSequence(seq, 1, ref, 1)
end

function AlignedSequence(seq::Sequence, seqpos::Integer,
                         ref::Sequence, refpos::Integer)
    if length(seq) != length(ref)
        throw(ArgumentError("two sequences must be the same length"))
    end
    seqpos -= 1
    refpos -= 1
    op = OP_START
    newseq = similar(seq, 0)  # sequence without gap symbols
    anchors = AlignmentAnchor[]
    for (x, y) in zip(seq, ref)
        if x == gap(eltype(seq)) && y == gap(eltype(ref))
            throw(ArgumentError("two sequences must not have gaps at the same position"))
        elseif x == gap(eltype(seq))
            op′ = OP_DELETE
        elseif y == gap(eltype(ref))
            op′ = OP_INSERT
        elseif x == y
            op′ = OP_SEQ_MATCH
        else
            op′ = OP_SEQ_MISMATCH
        end

        if op′ != op
            push!(anchors, AlignmentAnchor(seqpos, refpos, op))
            op = op′
        end

        if x != gap(eltype(seq))
            seqpos += 1
            push!(newseq, x)
        end
        if y != gap(eltype(ref))
            refpos += 1
        end
    end
    push!(anchors, AlignmentAnchor(seqpos, refpos, op))
    return AlignedSequence(newseq, anchors)
end

"""
First position in the reference sequence.
"""
function IntervalTrees.first(alnseq::AlignedSequence)
    return alnseq.aln.firstref
end

"""
Last position in the reference sequence.
"""
function IntervalTrees.last(alnseq::AlignedSequence)
    return alnseq.aln.lastref
end

function seq2ref(i::Integer, alnseq::AlignedSequence)
    return seq2ref(i, alnseq.aln)
end

function ref2seq(i::Integer, alnseq::AlignedSequence)
    return ref2seq(i, alnseq.aln)
end

# simple letters and dashes representation of an alignment
function Base.show(io::IO, alnseq::AlignedSequence)
    # print a representation of the reference sequence
    anchors = alnseq.aln.anchors
    for i in 2:length(anchors)
        if ismatchop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '·')
            end
        elseif isinsertop(anchors[i].op)
            for _ in anchors[i-1].seqpos+1:anchors[i].seqpos
                write(io, '-')
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '·')
            end
        end
    end
    write(io, '\n')

    for i in 2:length(anchors)
        if ismatchop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                print(io, alnseq.seq[i])
            end
        elseif isinsertop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                print(io, alnseq.seq[i])
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '-')
            end
        end
    end
end
