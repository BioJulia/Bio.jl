# conversion between different alphabet size
for A in [DNAAlphabet, RNAAlphabet]
    # 4 bits => 2 bits
    @eval function Base.convert(::Type{BioSequence{$(A{2})}}, seq::BioSequence{$(A{4})})
        # TODO: make it faster with bit-parallel algorithm
        newseq = BioSequence{$(A{2})}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end

    # 2 bits => 4 bits
    @eval function Base.convert(::Type{BioSequence{$(A{4})}}, seq::BioSequence{$(A{2})})
        newseq = BioSequence{$(A{4})}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end
end

# conversion between DNA and RNA
for (A1, A2) in [(DNAAlphabet, RNAAlphabet), (RNAAlphabet, DNAAlphabet)], n in (2, 4)
    # NOTE: assumes that binary representation is identical between DNA and RNA
    @eval function Base.convert(::Type{BioSequence{$(A1{n})}},
                                seq::BioSequence{$(A2{n})})
        newseq = BioSequence{$(A1{n})}(seq.data, seq.part, true)
        seq.shared = true
        return newseq
    end
end

# from a vector
function Base.convert{A<:DNAAlphabet}(::Type{BioSequence{A}},
                                      seq::AbstractVector{DNA})
    return BioSequence{A}(seq, 1, endof(seq))
end
function Base.convert{A<:RNAAlphabet}(::Type{BioSequence{A}},
                                      seq::AbstractVector{RNA})
    return BioSequence{A}(seq, 1, endof(seq))
end
function Base.convert(::Type{AminoAcidSequence}, seq::AbstractVector{AminoAcid})
    return AminoAcidSequence(seq, 1, endof(seq))
end

# to a vector
Base.convert(::Type{Vector}, seq::BioSequence) = collect(seq)
function Base.convert{A<:DNAAlphabet}(::Type{Vector{DNA}},
                                      seq::BioSequence{A})
    return collect(seq)
end
function Base.convert{A<:RNAAlphabet}(::Type{Vector{RNA}},
                                      seq::BioSequence{A})
    return collect(seq)
end
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = collect(seq)

# from/to a string
function Base.convert{S<:AbstractString}(::Type{S}, seq::BioSequence)
    return convert(S, [Char(x) for x in seq])
end
Base.convert{S<:AbstractString,A}(::Type{BioSequence{A}}, seq::S) = BioSequence{A}(seq)
