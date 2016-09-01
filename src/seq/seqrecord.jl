# Sequence Record
# ===============
#
# Sequence record with name and attached metadata.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
`SeqRecord{S,T}` is a type holding a named sequence of type `S`, along with
some arbitrary metadata of type `T`.
"""
type SeqRecord{S<:Sequence,T}
    name::StringField
    seq::S
    metadata::T
end

function SeqRecord(name::AbstractString, seq::Sequence, metadata=nothing)
    return SeqRecord(StringField(name), seq, metadata)
end

function (::Type{SeqRecord{S,T}}){S,T}()
    return SeqRecord{S,T}(StringField(), S(), T())
end

"""
    seqname(record)

Return the sequence name of `record`.
"""
function seqname(rec::SeqRecord)
    return rec.name
end

"""
    sequence(record)

Return the sequence of `record`.
"""
function sequence(rec::SeqRecord)
    return rec.seq
end

"""
    metadata(record)

Return the metadata of `record`.
"""
function metadata(rec::SeqRecord)
    return rec.metadata
end

function seqtype{S,T}(::Type{SeqRecord{S,T}})
    return S
end

function Base.:(==){T<:SeqRecord}(a::T, b::T)
    return a.name == b.name && a.seq == b.seq && a.metadata == b.metadata
end

function Base.copy{T <: SeqRecord}(seqrec::T)
    if seqrec.metadata == nothing
        # no copy method for the Void type
        return T(copy(seqrec.name), copy(seqrec.seq), nothing)
    else
        return T(copy(seqrec.name), copy(seqrec.seq), copy(seqrec.metadata))
    end
end

function Base.show(io::IO, rec::SeqRecord)
    println(io, summary(rec), ':')
    println(io, "  name: ", seqname(rec))
    println(io, "  sequence: ", string_compact(sequence(rec)))
      print(io, "  metadata: ", metadata(rec))
end
