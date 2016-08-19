# SAM Record
# ==========

type SAMRecord
    name::Bio.StringFields.StringField
    flag::UInt16
    refname::Bio.StringFields.StringField
    pos::Int64
    mapq::Int8
    cigar::Bio.StringFields.StringField
    next_refname::Bio.StringFields.StringField
    next_pos::Int64
    tlen::Int32
    seq::DNASequence
    qual::Vector{UInt8}
    optional_fields::Dict{String,Any}
end

function SAMRecord()
    return SAMRecord("", 0x0000, "*", 0, 0, "*", "*", 0, 0, "", UInt8[], Dict())
end

function Base.isless(rec1::SAMRecord, rec2::SAMRecord)
    # left-most position of alignment
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

function Base.show(io::IO, rec::SAMRecord)
    println(summary(rec), ':')
    println(io, "reference name: ", refname(rec))
    println(io, "next reference name: ", nextrefname(rec))
    println(io, "position: ", position(rec))
    println(io, "next position: ", nextposition(rec))
    println(io, "mapping quality: ", mappingquality(rec))
    println(io, "flag: ", flag(rec))
    println(io, "template length: ", templatelength(rec))
    println(io, "sequence name: ", seqname(rec))
    println(io, "CIGAR string: ", cigar(rec))
    println(io, "sequence: ", sequence(rec))
    println(io, "base qualities: ", qualities(rec))
      print(io, "optional fields: ", rec.optional_fields)
end

function Base.copy(rec::SAMRecord)
    return deepcopy(rec)
end

function ismapped(r::SAMRecord)
    return r.pos != 0
end

function Bio.Intervals.seqname(r::SAMRecord)
    return r.name
end

function flag(r::SAMRecord)
    return r.flag
end

function refname(r::SAMRecord)
    return r.refname
end

function nextrefname(r::SAMRecord)
    return r.next_refname
end

function Base.position(r::SAMRecord)
    return r.pos
end

function nextposition(r::SAMRecord)
    return r.next_pos
end

function mappingquality(r::SAMRecord)
    return r.mapq
end

function templatelength(r::SAMRecord)
    return r.tlen
end

function cigar(r::SAMRecord)
    return r.cigar
end

function sequence(r::SAMRecord)
    return r.seq
end

function qualities(r::SAMRecord)
    return r.qual
end

function Base.getindex(r::SAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return r.optional_fields[tag]
end

function Base.setindex!(r::SAMRecord, val, tag::AbstractString)
    checkkeytag(tag)
    setindex!(r.optional_fields, val, tag)
    return r
end

function Base.delete!(r::SAMRecord, tag::AbstractString)
    checkkeytag(tag)
    delete!(r.optional_fields, tag)
    return r
end

function Base.haskey(r::SAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return haskey(r.optional_fields, tag)
end
