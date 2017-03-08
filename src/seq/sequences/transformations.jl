# Transformations
# ---------------

function Base.push!{A}(seq::BioSequence{A}, x)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    encoded_setindex!(seq, bin, endof(seq))
    return seq
end

function Base.pop!(seq::BioSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    x = seq[end]
    deleteat!(seq, endof(seq))
    return x
end

function Base.insert!{A}(seq::BioSequence{A}, i::Integer, x)
    checkbounds(seq, i)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    copy!(seq, i + 1, seq, i, endof(seq) - i)
    encoded_setindex!(seq, bin, i)
    return seq
end

function Base.deleteat!{A,T<:Integer}(seq::BioSequence{A}, range::UnitRange{T})
    checkbounds(seq, range)
    copy!(seq, range.start, seq, range.stop + 1, length(seq) - range.stop)
    resize!(seq, length(seq) - length(range))
    return seq
end

function Base.deleteat!(seq::BioSequence, i::Integer)
    checkbounds(seq, i)
    copy!(seq, i, seq, i + 1, length(seq) - i)
    resize!(seq, length(seq) - 1)
    return seq
end

function Base.append!{A}(seq::BioSequence{A}, other::BioSequence{A})
    resize!(seq, length(seq) + length(other))
    copy!(seq, endof(seq) - length(other) + 1, other, 1)
    return seq
end

function Base.shift!(seq::BioSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    x = seq[1]
    deleteat!(seq, 1)
    return x
end

function Base.unshift!{A}(seq::BioSequence{A}, x)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    copy!(seq, 2, seq, 1, length(seq) - 1)
    encoded_setindex!(seq, bin, 1)
    return seq
end

"""
    resize!(seq, size)

Resize a biological sequence `seq`, to a given `size`.
"""
function Base.resize!{A}(seq::BioSequence{A}, size::Integer)
    if size < 0
        throw(ArgumentError("size must be non-negative"))
    end
    orphan!(seq, size)
    resize!(seq.data, seq_data_len(A, size + seq.part.start - 1))
    seq.part = seq.part.start:seq.part.start+size-1
    return seq
end

"""
    empty!(seq)

Completely empty a biological sequence `seq` of nucleotides.
"""
Base.empty!(seq::BioSequence) = resize!(seq, 0)

function Base.filter!{A}(f::Function, seq::BioSequence{A})
    orphan!(seq)

    len = 0
    next = bitindex(seq, 1)
    j = index(next)
    datum::UInt64 = 0
    for i in 1:endof(seq)
        x = inbounds_getindex(seq, i)
        if f(x)
            datum |= enc64(seq, x) << offset(next)
            len += 1
            next += bitsof(A)
            if index(next) != j
                seq.data[j] = datum
                datum = 0
                j = index(next)
            end
        end
    end
    if offset(next) > 0
        seq.data[j] = datum
    end
    resize!(seq, len)

    return seq
end

function Base.filter(f::Function, seq::BioSequence)
    return filter!(f, copy(seq))
end

function Base.map!(f::Function, seq::BioSequence)
    orphan!(seq)
    for i in 1:endof(seq)
        unsafe_setindex!(seq, f(inbounds_getindex(seq, i)), i)
    end
    return seq
end

function Base.map(f::Function, seq::BioSequence)
    return map!(f, copy(seq))
end

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

# Shuffle
# -------

function Base.shuffle(seq::BioSequence)
    return shuffle!(copy(seq))
end

function Base.shuffle!(seq::BioSequence)
    orphan!(seq)
    # Fisher-Yates shuffle
    for i in 1:endof(seq)-1
        j = rand(i:endof(seq))
        seq[i], seq[j] = seq[j], seq[i]
    end
    return seq
end
