type AlignedSequence
  Name::String
  GapSpace::Vector{Int}
  Mismatches::Vector{Int}
  MismatchSpace::Vector{Int}
end

type Alignment{T}
  Reference::{T}
  AlignedSequences::Vector{AlignedSequence}
end