# BioSequence
# ===========
#
# A general purpose biological sequence representation.
#
# It is cheap to create a subsequence from a sequence because sequences can
# share underlying data: creating a subsequence is just allocating a new
# sequence object defining a part of the original sequence without copy.
# Destructive operations create a copy of the underlying data if and only if it
# is shared between (sub)sequences. This is often called as copy-on-write
# strategy in computer science and should be transparent to the user.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

# About Internals
# ---------------
#
# The `data` field of a `BioSequence{A}` object contains binary representation
# of a biological character sequence. Each character is encoded with an encoder
# corresponding to the alphabet `A` and compactly packed into `data`. To extract
# a character from a sequence, you should decode this binary sequence with a
# decoder that is a pair of the encoder. The length of encoded binary bits is
# fixed, and hence a character at arbitrary position can be extracted in a
# constant time. To know the exact location of a character at a position, you
# can use the `bitsof(seq, i)` function, which returns a pair of element's index
# containing binary bits and bits' offset. As a whole, character extraction
# `seq[i]` can be written as:
#
#     j = bitindex(seq, i)
#     decode(A, (seq.data[index(j)] >> offset(j)) & mask(A))
#
#  index :        index(j) - 1       index(j)       index(j) + 1
#   data : ....|xxxxx...........|xxXxxxxxxxxxxxxx|............xxxx|....
# offset :                          |<-offset(j)-|
#  width :      |<---- 64 ---->| |<---- 64 ---->| |<---- 64 ---->|
#
#  * '.' : unused (4 bits/char)
#  * 'x' : used
#  * 'X' : used and pointed by index `i`

"""
Biological sequence data structure indexed by an alphabet type `A`.
"""
type BioSequence{A<:Alphabet} <: Sequence
    data::Vector{UInt64}  # encoded character sequence data
    part::UnitRange{Int}  # interval within `data` defining the (sub)sequence
    shared::Bool          # true if and only if `data` is shared between sequences

    function BioSequence(data::Vector{UInt64}, part::UnitRange{Int}, shared::Bool)
        return new(data, part, shared)
    end
end

typealias DNASequence       BioSequence{DNAAlphabet{4}}
typealias RNASequence       BioSequence{RNAAlphabet{4}}
typealias AminoAcidSequence BioSequence{AminoAcidAlphabet}
typealias CharSequence      BioSequence{CharAlphabet}


function seq_data_len{A}(::Type{A}, len::Integer)
    return cld(len, div(64, bitsof(A)))
end


# Summaries
# ---------

Base.summary{A<:DNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "DNA Sequence")
Base.summary{A<:RNAAlphabet}(seq::BioSequence{A}) = string(length(seq), "nt ", "RNA Sequence")
Base.summary(seq::AminoAcidSequence) = string(length(seq), "aa ", "Amino Acid Sequence")
Base.summary(seq::CharSequence) = string(length(seq), "char ", "Char Sequence")
