# MurmurHash3
# ===========
#
# The hash function defined here cares about the starting position of a
# character sequence in the underlying data. That means, even if the starting
# positions of two sequences (`s1` and `s2`) are different in their `data`
# field, their hash values are identical if `s1 == s1` is true.
#
# The implementation is based on the 128bit MurmurHash3 function, which was
# written by Austin Appleby, and the source code is distributed under the public
# domain: https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp

@inline function rotl64(x::UInt64, r)
    return (x << r) | (x >> (64 - r))
end

@inline function fmix64(k::UInt64)
    k $= k >> 33
    k *= 0xff51afd7ed558ccd
    k $= k >> 33
    k *= 0xc4ceb9fe1a85ec53
    k $= k >> 33
    return k
end

macro murmur1()
    esc(quote
        k1 *= c1
        k1 = rotl64(k1, 31)
        k1 *= c2
        h1 $= k1
    end)
end

macro murmur2()
    esc(quote
        k2 *= c2
        k2 = rotl64(k2, 33)
        k2 *= c1
        h2 $= k2
    end)
end

macro murmur()
    esc(quote
        @murmur1
        h1 = rotl64(h1, 27)
        h1 += h2
        h1 = h1 * 5 + 0x52dce729

        @murmur2
        h2 = rotl64(h2, 31)
        h2 += h1
        h2 = h2 * 5 + 0x38495ab5
    end)
end

# ref: MurmurHash3_x64_128
function Base.hash(seq::BioSequence, seed::UInt64)
    # Mix sequence length so that dna"A" and dna"AA"
    # return the different hash values.
    h1::UInt64 = h2::UInt64 = hash(length(seq), seed)
    c1 = 0x87c37b91114253d5
    c2 = 0x4cf5ad432745937f

    next = bitindex(seq, 1)
    last = bitindex(seq, endof(seq) + 1)

    k1::UInt64 = 0
    k2::UInt64 = 0

    # body
    r = offset(next)
    data = seq.data
    if last - next ≥ 128
        if r == 0
            @inbounds while last - next ≥ 128
                j = index(next)
                k1 = data[j]
                k2 = data[j+1]
                @murmur
                next += 128
            end
        else
            x = data[index(next)]
            @inbounds while last - next ≥ 128
                j = index(next)
                y = data[j+1]
                z = data[j+2]
                k1 = x >> r | y << (64 - r)
                k2 = y >> r | z << (64 - r)
                @murmur
                x = z
                next += 128
            end
        end
    end

    # tail
    k1 = 0
    k2 = 0
    if next < last
        x = data[index(next)]
        k1 |= x >> r
        m1 = mask(last - next)
        m2 = mask(max(last - (next + 64), 0))
        next += 64 - r
        if next < last
            y = data[index(next)]
            k1 |= y << (64 - r)
            k2 |= y >> r
            next += 64
            if next < last
                z = data[index(next)]
                k2 |= z << (64 - r)
            end
        end
        k1 &= m1
        k2 &= m2
        @murmur1
        @murmur2
    end

    # finalization
    h1 $= length(seq)
    h2 $= length(seq)
    h1 += h2
    h2 += h1
    h1 = fmix64(h1)
    h2 = fmix64(h2)
    h1 += h2
    h2 += h1

    return h1
end
