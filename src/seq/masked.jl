# Masked Sequence
# ===============

"""
Masked sequence data structure.
"""
immutable MaskedSequence{S<:Sequence,M<:AbstractVector{Bool}} <: Sequence
    seq::S
    mask::M

    function MaskedSequence(seq::Sequence, mask::AbstractVector{Bool})
        if length(seq) != length(mask)
            throw(ArgumentError("mask length doesn't match"))
        end
        # create a subsequence
        seq = seq[1:end]
        return new(seq, mask)
    end
end

function MaskedSequence(seq::Sequence, mask::AbstractVector{Bool})
    return MaskedSequence{typeof(seq),typeof(mask)}(seq, mask)
end

function mask(p::Function, seq::Sequence)
    mask = falses(length(seq))
    for (i, x) in enumerate(seq)
        if p(x)
            mask[i] = true
        end
    end
    return MaskedSequence(seq, mask)
end

function ismasked(seq::MaskedSequence, i::Integer)
    return seq.mask[i]
end

function Base.length(seq::MaskedSequence)
    return length(seq.seq)
end

function Base.eltype(seq::MaskedSequence)
    return eltype(seq.seq)
end

function inbounds_getindex(seq::MaskedSequence, i::Integer)
    return inbounds_getindex(seq.seq, i)
end
