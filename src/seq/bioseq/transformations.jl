# Transformations
# ---------------

"""
    reverse!(seq)

Reverse a biological sequence `seq` in place.
"""
function Base.reverse!(seq::BioSequence)
    orphan!(seq)
    for i in 1:div(endof(seq), 2)
        x = inbounds_getindex(seq, i)
        unsafe_setindex!(seq, inbounds_getindex(seq, endof(seq) - i + 1), i)
        unsafe_setindex!(seq, x, endof(seq) - i + 1)
    end
    return seq
end

"""
    reverse(seq)

Create a sequence which is the reverse of the bioloigcal sequence `seq`.
"""
Base.reverse(seq::BioSequence) = reverse!(copy(seq))

@generated function Base.reverse{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    n = bitsof(A)
    if n == 2
        nucrev = :nucrev2
    elseif n == 4
        nucrev = :nucrev4
    else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        data = Vector{UInt64}(seq_data_len(A, length(seq)))
        i = 1
        next = bitindex(seq, endof(seq))
        stop = bitindex(seq, 0)
        r = rem(offset(next) + $n, 64)
        if r == 0
            @inbounds while next - stop > 0
                x = seq.data[index(next)]
                data[i] = $nucrev(x)
                i += 1
                next -= 64
            end
        else
            @inbounds while next - stop > 64
                j = index(next)
                x = (seq.data[j] << (64 - r)) | (seq.data[j-1] >> r)
                data[i] = $nucrev(x)
                i += 1
                next -= 64
            end
            if next - stop > 0
                j = index(next)
                x = seq.data[j] << (64 - r)
                if r < next - stop
                    x |= seq.data[j-1] >> r
                end
                data[i] = $nucrev(x)
            end
        end
        return BioSequence{A}(data, 1:length(seq), false)
    end
end

@inline function nucrev2(x::UInt64)
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>  2
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
     return x
end

@inline function nucrev4(x::UInt64)
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
     return x
end

"""
    complement!(seq)

Make a complement sequence of `seq` in place.
"""
function complement!{A<:Union{DNAAlphabet{2},RNAAlphabet{2}}}(seq::BioSequence{A})
    orphan!(seq)
    next = bitindex(seq, 1)
    stop = bitindex(seq, endof(seq) + 1)
    @inbounds while next < stop
        seq.data[index(next)] = ~seq.data[index(next)]
        next += 64
    end
    return seq
end

"""
    complement!(seq)

Transform `seq` into it's complement.
"""
function complement!{A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}(seq::BioSequence{A})
    orphan!(seq)
    next = bitindex(seq, 1)
    stop = bitindex(seq, endof(seq) + 1)
    @inbounds while next < stop
        x = seq.data[index(next)]
        seq.data[index(next)] = (
            ((x & 0x1111111111111111) << 3) | ((x & 0x8888888888888888) >> 3) |
            ((x & 0x2222222222222222) << 1) | ((x & 0x4444444444444444) >> 1))
        next += 64
    end
    return seq
end

"""
    complement(seq)

Make a complement sequence of `seq`.
"""
function complement{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(copy(seq))
end

"""
    reverse_complement!(seq)

Make a reversed complement sequence of `seq` in place.

Ambiguous nucleotides are left as-is.
"""
function reverse_complement!{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(reverse!(seq))
end

"""
    reverse_complement(seq)

Make a reversed complement sequence of `seq`.

Ambiguous nucleotides are left as-is.
"""
function reverse_complement{A<:Union{DNAAlphabet,RNAAlphabet}}(seq::BioSequence{A})
    return complement!(reverse(seq))
end
