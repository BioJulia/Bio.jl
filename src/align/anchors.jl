# Alignment Anchor
# ================
#
# Sequence alignment anchor type and alignment data data structures.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
A type to store the operation enocded in an alignment (CIGAR Operations).
Also stores the position in the alignment view of the sequences, and the
corresponding position in the unaltered source sequence (nucleotide or protein).

It stores the position values as Int types and the alignment operation is stored
as a type `Operation`, these are defined in the file `operations.jl` and loaded
into Bio.Align namespace as a series of global constants.
"""
immutable AlignmentAnchor
    seqpos::Int
    refpos::Int
    op::Operation
end

function AlignmentAnchor(pos::Tuple{Int,Int}, op)
    return AlignmentAnchor(pos[1], pos[2], op)
end

function Base.show(io::IO, anc::AlignmentAnchor)
    print(io, "AlignmentAnchor(", anc.seqpos, ", ", anc.refpos, ", '", anc.op, "')")
end


# Alignment
# ---------

immutable Alignment
    anchors::Vector{AlignmentAnchor}
    firstref::Int
    lastref::Int

    function Alignment(anchors::Vector{AlignmentAnchor}, check::Bool=true)
        # optionally check coherence of the anchors
        if check
            # empty alignments are valid, representing an unaligned sequence
            if !isempty(anchors)
                if anchors[1].op != OP_START
                    error("Alignments must begin with on OP_START anchor.")
                end

                for i in 2:length(anchors)
                    if anchors[i].refpos < anchors[i-1].refpos ||
                       anchors[i].seqpos < anchors[i-1].seqpos
                        error("Alignment anchors must be sorted.")
                    end

                    op = anchors[i].op
                    if convert(UInt8, op) > convert(UInt8, OP_MAX_VALID)
                        error("Anchor at index $(i) has an invalid operation.")
                    end

                    # reference skip/delete operations
                    if isdeleteop(op)
                        if anchors[i].seqpos != anchors[i-1].seqpos
                            error("Invalid anchor positions for reference deletion.")
                        end
                    # reference insertion operations
                    elseif isinsertop(op)
                        if anchors[i].refpos != anchors[i-1].refpos
                            error("Invalid anchor positions for reference insertion.")
                        end
                    # match operations
                    elseif ismatchop(op)
                        if anchors[i].refpos - anchors[i-1].refpos !=
                           anchors[i].seqpos - anchors[i-1].seqpos
                            error("Invalid anchor positions for match operation.")
                        end
                    end
                end
            end
        end

        # compute first and last aligned reference positions
        firstref = 0
        for i in 1:length(anchors)
            if ismatchop(anchors[i].op)
                firstref = anchors[i-1].refpos + 1
                break
            end
        end

        lastref = 0
        for i in length(anchors):-1:1
            if ismatchop(anchors[i].op)
                lastref = anchors[i].refpos
                break
            end
        end

        return new(anchors, firstref, lastref)
    end
end

"""
Construct an `Alignment`
"""
function Alignment(cigar::AbstractString, seqpos::Int=1, refpos::Int=1)
    # path starts prior to the first aligned position pair
    seqpos -= 1
    refpos -= 1

    n = 0
    anchors = AlignmentAnchor[AlignmentAnchor(seqpos, refpos, OP_START)]
    for c in cigar
        if isdigit(c)
            n = n * 10 + convert(Int, c - '0')
        else
            if n == 0
                error("CIGAR operations must be prefixed by a positive integer.")
            end
            op = Operation(c)
            if ismatchop(op)
                seqpos += n
                refpos += n
            elseif isinsertop(op)
                seqpos += n
            elseif isdeleteop(op)
                refpos += n
            else
                error("The $(op) CIGAR operation is not yet supported.")
            end

            push!(anchors, AlignmentAnchor(seqpos, refpos, op))
            n = 0
        end
    end

    return Alignment(anchors)
end

function Base.(:(==))(a::Alignment, b::Alignment)
    return a.anchors == b.anchors && a.firstref == b.firstref && a.lastref == b.lastref
end

function Base.show(io::IO, aln::Alignment)
    # print a representation of the reference sequence
    anchors = aln.anchors
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

    # print a representation of the aligned sequence
    for i in 2:length(anchors)
        if ismatchop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                write(io, '·')
            end
        elseif isinsertop(anchors[i].op)
            for i in anchors[i-1].seqpos+1:anchors[i].seqpos
                write(io, '·')
            end
        elseif isdeleteop(anchors[i].op)
            for _ in anchors[i-1].refpos+1:anchors[i].refpos
                write(io, '-')
            end
        end
    end
end

# find the index of the first anchor that satisfies `i ≤ pos`
@generated function findanchor{isseq}(aln::Alignment, i::Integer, ::Type{Val{isseq}})
    pos = isseq ? :seqpos : :refpos
    quote
        anchors = aln.anchors
        lo = 1
        hi = endof(anchors)
        if !(anchors[lo].$pos < i ≤ anchors[hi].$pos)
            return 0
        end
        # binary search
        @inbounds while hi - lo > 2
            m = (lo + hi) >> 1
            if anchors[m].$pos < i
                lo = m
            else  # i ≤ anchors[m].$pos
                hi = m
            end
            # invariant (activate this for debugging)
            #@assert anchors[lo].$pos < i ≤ anchors[hi].$pos
        end
        # linear search
        @inbounds for j in lo+1:hi
            if i ≤ aln.anchors[j].$pos
                return j
            end
        end
        # do not reach here
        @assert false
        return 0
    end
end

"""
    seq2ref(i, aln)

Map a position from sequence to reference.
"""
function seq2ref(i::Integer, aln::Alignment)
    idx = findanchor(aln, i, Val{true})
    if idx == 0
        throw(ArgumentError("invalid sequence position: $i"))
    end
    anchor = aln.anchors[idx]
    refpos = anchor.refpos
    if ismatchop(anchor.op)
        refpos += i - anchor.seqpos
    end
    return refpos, anchor.op
end

"""
    ref2seq(i, aln)

Map a position from reference to sequence.
"""
function ref2seq(i::Integer, aln::Alignment)
    idx = findanchor(aln, i, Val{false})
    if idx == 0
        throw(ArgumentError("invalid reference position: $i"))
    end
    anchor = aln.anchors[idx]
    seqpos = anchor.seqpos
    if ismatchop(anchor.op)
        seqpos += i - anchor.refpos
    end
    return seqpos, anchor.op
end


# AlignedSequence
# ---------------

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

"""
Output a CIGAR string encoding of an `Alignment`. This is not entirely lossless as
it discards the alignments start positions.
"""
function cigar(aln::Alignment)
    anchors = aln.anchors
    out = IOBuffer()
    seqpos = anchors[1].seqpos
    refpos = anchors[1].refpos

    for i in 2:length(anchors)
        n = max(anchors[i].seqpos - anchors[i-1].seqpos,
                anchors[i].refpos - anchors[i-1].refpos)
        print(out, n, Char(anchors[i].op))
    end
    return takebuf_string(out)
end

