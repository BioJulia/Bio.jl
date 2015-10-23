
# Anchor definition
# ------------------
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


# Basic operators for AlignmentAnchors
# -------------------------------------

function show(io::IO, anc::AlignmentAnchor)
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


function ==(a::Alignment, b::Alignment)
    return a.anchors == b.anchors && a.firstref == b.firstref && a.lastref == b.lastref
end


function show(io::IO, aln::Alignment)
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


immutable AlignedSequence{S}
    seq::S
    aln::Alignment
end


function AlignedSequence{S}(seq::S, aln::Alignment)
    return AlignedSequence{S}(seq, aln)
end


function AlignedSequence{S}(seq::S, anchors::Vector{AlignmentAnchor},
                            check::Bool=true)
    return AlignedSequence(seq, Alignment(anchors, check))
end


"""
First position in the reference sequence.
"""
function first(alnseq::AlignedSequence)
    return alnseq.aln.firstref
end


"""
Last position in the reference sequence.
"""
function last(alnseq::AlignedSequence)
    return alnseq.aln.lastref
end


# simple letters and dashes representation of an alignment
function show(io::IO, alnseq::AlignedSequence)
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

