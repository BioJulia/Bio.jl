# BitIndex
# --------
#
# Utils for indexing bits in a vector of 64-bit integers (internal use only).
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

immutable BitIndex
    val::Int64
end

BitIndex(index, nbits) = BitIndex((index - 1) << trailing_zeros(nbits))

#         index(i)-1        index(i)        index(i)+1
# ....|................|..X.............|................|....
#                          |<-offset(i)-|
#                      |<--- 64 bits -->|

index(i::BitIndex) = (i.val >> 6) + 1
offset(i::BitIndex) = i.val & 0b111111

@compat begin
    Base.:+(i::BitIndex, n::Int) = BitIndex(i.val + n)
    Base.:-(i::BitIndex, n::Int) = BitIndex(i.val - n)
    Base.:-(i1::BitIndex, i2::BitIndex) = i1.val - i2.val
    Base.:(==)(i1::BitIndex, i2::BitIndex) = i1.val == i2.val
end
Base.isless(i1::BitIndex, i2::BitIndex) = isless(i1.val, i2.val)
Base.cmp(i1::BitIndex, i2::BitIndex) = cmp(i1.val, i2.val)
Base.start(i::BitIndex) = 1
Base.done(i::BitIndex, s) = s > 2
Base.next(i::BitIndex, s) = ifelse(s == 1, (index(i), 2), (offset(i), 3))
Base.show(io::IO, i::BitIndex) = print(io, '(', index(i), ", ", offset(i), ')')
