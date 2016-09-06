# Alignment Anchor
# ================
#
# Sequence alignment anchor type.
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
