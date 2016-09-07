# Alignment
# =========
#
# Sequence alignment type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable Alignment
    anchors::Vector{AlignmentAnchor}
    firstref::Int
    lastref::Int

    function Alignment(anchors::Vector{AlignmentAnchor}, check::Bool=true)
        # optionally check coherence of the anchors
        if check
            check_alignment_anchors(anchors)
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

# Make an alignment object from a CIGAR string.
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

function Base.:(==)(a::Alignment, b::Alignment)
    return a.anchors == b.anchors && a.firstref == b.firstref && a.lastref == b.lastref
end

function Base.show(io::IO, aln::Alignment)
    println(io, summary(aln), ':')
    println(io, "  alignment: ", cigar(aln))
      print(io, "  aligned range: ", aln.firstref, '-', aln.lastref)
end

"""
    seq2ref(aln, i)

Map a position from sequence to reference.
"""
function seq2ref(aln::Alignment, i::Integer)
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
    ref2seq(aln, i)

Map a position from reference to sequence.
"""
function ref2seq(aln::Alignment, i::Integer)
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

"""
    cigar(aln::Alignment)

Make a CIGAR string encoding of `aln`.

This is not entirely lossless as it discards the alignments start positions.
"""
function cigar(aln::Alignment)
    anchors = aln.anchors
    if isempty(anchors)
        return ""
    end
    seqpos = anchors[1].seqpos
    refpos = anchors[1].refpos
    @assert anchors[1].op == OP_START
    out = IOBuffer()
    for i in 2:length(anchors)
        n = max(anchors[i].seqpos - anchors[i-1].seqpos,
                anchors[i].refpos - anchors[i-1].refpos)
        print(out, n, Char(anchors[i].op))
    end
    return takebuf_string(out)
end

# Check validity of a sequence of anchors.
function check_alignment_anchors(anchors)
    if isempty(anchors)
        # empty alignments are valid, representing an unaligned sequence
        return
    end

    if anchors[1].op != OP_START
        error("Alignments must begin with on OP_START anchor.")
    end

    for i in 2:endof(anchors)
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
