srcSeq = ['A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C']
gapanchors = [
    AlignmentAnchor(gp::Int = 7, sp::Int = 4, op::Operation = OP_GAP),
    AlignmentAnchor(gp::Int = 14, sp::Int = 9, op::Operation = OP_GAP)
    ]

matchAndMismatch = [
    AlignmentAnchor(gp::Int = 3, sp::Int = 3, op::Operation = OP_CIGAR_M),
    AlignmentAnchor(gp::Int = 7, sp::Int = 4, op::Operation = OP_CIGAR_N),
    AlignmentAnchor(gp::Int = 11, sp::Int = 8, op::Operation = OP_CIGAR_M),
    AlignmentAnchor(gp::Int = 14, sp::Int = 9, op::Operation = OP_CIGAR_N),
    AlignmentAnchor(gp::Int = 15, sp::Int = 11, op::Operation = OP_CIGAR_M),
]
