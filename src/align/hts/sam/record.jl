# SAM Record
# ==========

type SAMRecord
    name::StringField
    flag::UInt16
    refname::StringField
    pos::Int64
    mapq::UInt8
    cigar::StringField
    next_refname::StringField
    next_pos::Int64
    tlen::Int32
    seq::StringField
    qual::StringField
    optional_fields::Dict{String,Any}
end

function SAMRecord()
    return SAMRecord("*", SAM_FLAG_UNMAP, "*", 0, 0, "*", "*", 0, 0, "*", "*", Dict())
end

function Base.isless(rec1::SAMRecord, rec2::SAMRecord)
    # compared by the left-most position of an alignment
    if rec1.name == rec2.name
        return isless(rec1.pos, rec2.pos)
    else
        return isless(rec1.name, rec2.name)
    end
end

function Base.:(==)(rec1::SAMRecord, rec2::SAMRecord)
    return (
        rec1.name            == rec2.name         &&
        rec1.flag            == rec2.flag         &&
        rec1.refname         == rec2.refname      &&
        rec1.pos             == rec2.pos          &&
        rec1.mapq            == rec2.mapq         &&
        rec1.cigar           == rec2.cigar        &&
        rec1.next_refname    == rec2.next_refname &&
        rec1.next_pos        == rec2.next_pos     &&
        rec1.tlen            == rec2.tlen         &&
        rec1.seq             == rec2.seq          &&
        rec1.qual            == rec2.qual         &&
        rec1.optional_fields == rec2.optional_fields)
end

function Base.copy(rec::SAMRecord)
    return SAMRecord(
        copy(rec.name),
        rec.flag,
        copy(rec.refname),
        rec.pos,
        rec.mapq,
        copy(rec.cigar),
        copy(rec.next_refname),
        rec.next_pos,
        rec.tlen,
        copy(rec.seq),
        copy(rec.qual),
        deepcopy(rec.optional_fields))
end

function ismapped(rec::SAMRecord)
    return flag(rec) & SAM_FLAG_UNMAP == 0
end

function Bio.Seq.seqname(rec::SAMRecord)
    return rec.name
end

function flag(rec::SAMRecord)
    return rec.flag
end

function refname(rec::SAMRecord)
    return rec.refname
end

function nextrefname(rec::SAMRecord)
    return rec.next_refname
end

function Bio.Intervals.leftposition(rec::SAMRecord)
    return rec.pos
end

function Bio.Intervals.rightposition(rec::SAMRecord)
    return leftposition(rec) + alignment_length(rec) - 1
end

function nextleftposition(rec::SAMRecord)
    return rec.next_pos
end

function mappingquality(rec::SAMRecord)
    return rec.mapq
end

function templatelength(rec::SAMRecord)
    return rec.tlen
end

function cigar(rec::SAMRecord)
    return rec.cigar
end

function alignment(rec::SAMRecord)
    if ismapped(rec)
        return Alignment(rec.cigar, 1, leftposition(rec))
    else
        return Alignment(AlignmentAnchor[])
    end
end

function Bio.Seq.sequence(rec::SAMRecord)
    return rec.seq
end

function seqlength(rec::SAMRecord)
    if rec.seq == "*"
        throw(ArgumentError("no sequence available"))
    end
    return length(rec.seq)
end

function qualities(rec::SAMRecord)
    return rec.qual
end

function Base.getindex(rec::SAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return rec.optional_fields[tag]
end

function Base.setindex!(rec::SAMRecord, val, tag::AbstractString)
    checkkeytag(tag)
    setindex!(rec.optional_fields, val, tag)
    return rec
end

function Base.delete!(rec::SAMRecord, tag::AbstractString)
    checkkeytag(tag)
    delete!(rec.optional_fields, tag)
    return rec
end

function Base.haskey(rec::SAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return haskey(rec.optional_fields, tag)
end

function optional_fields(rec::SAMRecord)
    return rec.optional_fields
end

# Return the length of alignment.
function alignment_length(rec::SAMRecord)
    if rec.cigar == "*"
        return 0
    end
    length = 0
    len = 0  # operation length
    for c in rec.cigar
        if isnumber(c)
            len = len * 10 + (c - '0')
        elseif isalpha(c)
            op = convert(Operation, c)
            if ismatchop(op) || isdeleteop(op)
                length += len
                len = 0
            end
        else
            error("invalid character in CIGAR: '$(c)'")
        end
    end
    return length
end
