## Conserved sites

"""
    create_nibble_mask(::Type{Conserved}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides are the same.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Conserved}, a::UInt64, b::UInt64)
    #certainmask = create_nibble_mask(Certain, a, b)
    return create_nibble_mask(a $ b, 0x0000000000000000) & certainmask
end


## Mutated sites

"""
    create_nibble_mask(::Type{Mutated}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides are different.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Mutated}, a::UInt64, b::UInt64)
    certainmask = create_nibble_mask(Certain, a, b)
    return (~create_nibble_mask(a $ b, 0x0000000000000000)) & certainmask
end

"""
    create_nibble_mask(::Type{Transition}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides constitute a transition mutation.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Transition}, a::UInt64, b::UInt64)
    ts_filter = (a $ b) & create_nibble_mask(Certain, a, b)
    return create_nibble_mask(ts_filter, 0xAAAAAAAAAAAAAAAA) |
           create_nibble_mask(ts_filter, 0x5555555555555555)
end

"""
    create_nibble_mask(::Type{Transversion}, a::UInt64, b::UInt64)

Create a mask of the nibbles in two chunks of
BioSequence{(DNA|RNA)Nucleotide{4}} data that represent sites where the encoded
nucleotides constitute a transition mutation.

**This is an internal method and should not be exported.**
"""
@inline function create_nibble_mask(::Type{Transversion}, a::UInt64, b::UInt64)
    diffs = (a $ b) & create_nibble_mask(Certain, a, b)
    return create_nibble_mask(diffs, 0x3333333333333333) |
           create_nibble_mask(diffs, 0x6666666666666666) |
           create_nibble_mask(diffs, 0xCCCCCCCCCCCCCCCC) |
           create_nibble_mask(diffs, 0x9999999999999999)

end






## Conserved sites

"""
    count_nibbles(::Type{Conserved}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
Conserved between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.
**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
Conserved. For example, 'A' and 'R', or 'A' and '-' will not be counted.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Conserved}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Conserved, a, b))
end


## Mutated sites

"""
    count_nibbles(::Type{Mutated}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
any mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}} data.
**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Mutated}, a::UInt64, b::UInt64)
    return count_one_nibbles(create_nibble_mask(Mutated, a, b))
end

"""
    count_nibbles(::Type{Transition}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
transition mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}}
data.
**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Transition}, a::UInt64, b::UInt64)
    pairdelmask = create_nibble_mask(Certain, a, b)
    a &= pairdelmask
    b &= pairdelmask
    diffs = a $ b
    ts_filtered = diffs $ 0xAAAAAAAAAAAAAAAA
    return count_one_nibbles(ts_filtered) + count_one_nibbles(~ts_filtered)
end

"""
    count_nibbles(::Type{Transversion}, a::UInt64, b::UInt64)
An _internal_ function, _not for export_, which will count the number of
transversion mutations between two chunks of BioSequence{(DNA|RNA)Nucleotide{4}}
data.
**Note:** Ambiguous cases or cases with gaps are ignored and not counted as
mutated. For example, 'A' and 'R', or 'A' and '-' will not be counted.
**This is an internal method and should not be exported.**
"""
@inline function count_nibbles(::Type{Transversion}, a::UInt64, b::UInt64)
    pairdelmask = create_nibble_mask(Certain, a, b)
    a &= pairdelmask
    b &= pairdelmask
    diffs = a $ b
    tv_filtered_a = diffs $ 0x3333333333333333
    tv_filtered_b = diffs $ 0x6666666666666666
    return count_one_nibbles(tv_filtered_a) + count_one_nibbles(~tv_filtered_a) +
           count_one_nibbles(tv_filtered_b) + count_one_nibbles(~tv_filtered_b)
end
