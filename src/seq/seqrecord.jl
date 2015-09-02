
# A sequence record is a named sequence with attached metadata.
"""
`SeqRecord{S, T}` is a type holding a named sequence of type `S`, along with
some arbitrary metadata of type `T`.
"""
type SeqRecord{S, T}
    name::StringField
    seq::S
    metadata::T

    function SeqRecord(name, seq::S, metadata::T)
        return new(name, seq, metadata)
    end

    function SeqRecord()
        return new("", S(), T())
    end
end


# Degelgate sequence operations
"Return a `SeqRecord` holding just the nucleotide at position `i`"
function Base.getindex(seqrec::SeqRecord, i::Integer)
    return SeqRecord(seqreq.name, seqrec.seq[i], seqreq.metadata)
end

"Return a `SeqRecord` holding the specified subsequence"
function Base.getindex(seqrec::SeqRecord, r::UnitRange)
    return SeqRecord(seqrec.name, seqrec.seq[r], seqrec.metadata)
end


typealias DNASeqRecord{T}       SeqRecord{DNASequence, T}
typealias RNASeqRecord{T}       SeqRecord{RNASequence, T}
typealias AminoAcidSeqRecord{T} SeqRecord{AminoAcidSequence, T}
