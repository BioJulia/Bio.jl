# BAM Record
# ==========

type BAMRecord
    # fixed-length fields (see BMA specs for the details)
    block_size::Int32
    refid::Int32
    pos::Int32
    bin_mq_nl::UInt32
    flag_nc::UInt32
    l_seq::Int32
    next_refid::Int32
    next_pos::Int32
    tlen::Int32

    # variable length data
    data::Vector{UInt8}

    # pointer to reference sequence names (shared)
    refseqnames::Vector{String}
end

# the data size of fixed-length fields (.block_size-.tlen)
const BAM_FIXED_FIELDS_BYTES = 36

function BAMRecord()
    flag_nc = UInt32(SAM_FLAG_UNMAP) << 16
    return BAMRecord(0, -1, -1, 0, flag_nc, 0, -1, -1, 0, UInt8[], String[])
end

# NOTE: this does not copy `refseqnames`.
function Base.copy(rec::BAMRecord)
    return BAMRecord(
        rec.block_size,
        rec.refid,
        rec.pos,
        rec.bin_mq_nl,
        rec.flag_nc,
        rec.l_seq,
        rec.next_refid,
        rec.next_pos,
        rec.tlen,
        copy(rec.data),
        rec.refseqnames)
end

function Base.isless(rec1::BAMRecord, rec2::BAMRecord)
    # compared by left-most position of an alignment
    if rec1.refid == rec2.refid
        return isless(rec1.pos, rec2.pos)
    else
        return isless(rec1.refid, rec2.refid)
    end
end

function Base.:(==)(rec1::BAMRecord, rec2::BAMRecord)
    return (
        rec1.refid      == rec2.refid      &&
        rec1.pos        == rec2.pos        &&
        rec1.bin_mq_nl  == rec2.bin_mq_nl  &&
        rec1.flag_nc    == rec2.flag_nc    &&
        rec1.l_seq      == rec2.l_seq      &&
        rec1.next_refid == rec2.next_refid &&
        rec1.next_pos   == rec2.next_pos   &&
        rec1.tlen       == rec2.tlen)  # TODO: check data
end

"""
    ismapped(rec::BAMRecord)

Return `true` if and only if `rec` is mapped to a reference sequence.
"""
function ismapped(rec::BAMRecord)
    return flag(rec) & SAM_FLAG_UNMAP == 0
end

"""
    refindex(rec::BAMRecord)

Return the index of a reference sequence that `rec` is mapped onto.

The index is 1-based and will be 0 for an alignment without mapping position.
"""
function refindex(rec::BAMRecord)
    return rec.refid + 1
end

function nextrefindex(rec::BAMRecord)
    return rec.next_refid + 1
end

"""
    refname(rec::BAMRecord)

Return the name of a reference sequence that `rec` is mapped onto.

If `rec` is unmapped, it returns `"*"` like SAM records.
"""
function refname(rec::BAMRecord)
    i = refindex(rec)
    if i == 0
        return "*"
    else
        return rec.refseqnames[i]
    end
end

function nextrefname(rec::BAMRecord)
    i = nextrefindex(rec)
    if i == 0
        return "*"
    else
        return rec.refseqnames[i]
    end
end

"""
    leftposition(rec::BAMRecord)

Return the leftmost mapping position of `rec`.

The index is 1-based and will be 0 for an alignment without mapping position.
"""
function Bio.Intervals.leftposition(rec::BAMRecord)
    return rec.pos + 1
end

"""
    rightposition(rec::BAMRecord)

Return the rightmost mapping position of `rec`.
"""
function Bio.Intervals.rightposition(rec::BAMRecord)
    return Int32(leftposition(rec) + alignment_length(rec) - 1)
end


function nextleftposition(rec::BAMRecord)
    return rec.next_pos + 1
end

"""
    mappingquality(rec::BAMRecord)

Return the mapping quality of the alignment `rec`.
"""
function mappingquality(rec::BAMRecord)
    return UInt8((rec.bin_mq_nl >> 8) & 0xff)
end

"""
    flag(rec::BAMRecord)

Return the flag of the alignment `rec`.
"""
function flag(rec::BAMRecord)
    return UInt16(rec.flag_nc >> 16)
end

"""
    templatelength(rec::BAMRecord)

Return the template length of the alignment `rec`.
"""
function templatelength(rec::BAMRecord)
    return rec.tlen
end

function Bio.Seq.seqname(rec::BAMRecord)
    # drop the last NUL character
    return unsafe_string(pointer(rec.data), max(seqname_length(rec) - 1, 0))
end

"""
    cigar_rle(rec::BAMRecord)

Return a run-length encoded tuple `(ops, lens)` of the CIGAR string.
See also `cigar`.
"""
function cigar_rle(rec::BAMRecord)
    offset = seqname_length(rec)
    ops = Bio.Align.Operation[]
    lens = Int[]
    for i in offset+1:4:offset+n_cigar_op(rec)*4
        x = unsafe_load(Ptr{UInt32}(pointer(rec.data, i)))
        op = Bio.Align.Operation(x & 0x0f)
        push!(ops, op)
        push!(lens, x >> 4)
    end
    return ops, lens
end

"""
    cigar(rec::BAMRecord)

Return a CIGAR string of the alignment `rec`. See also `cigar_rle`.
"""
function cigar(rec::BAMRecord)
    buf = IOBuffer()
    for (op, len) in zip(cigar_rle(rec)...)
        print(buf, len, Char(op))
    end
    return takebuf_string(buf)
end

"""
    alignment(rec::BAMRecord)

Make an alignment object from `rec`.
"""
function alignment(rec::BAMRecord)
    if !ismapped(rec)
        return Alignment(AlignmentAnchor[])
    end
    seqpos = 0
    refpos = leftposition(rec) - 1
    anchors = [AlignmentAnchor(seqpos, refpos, OP_START)]
    for (op, len) in zip(cigar_rle(rec)...)
        if ismatchop(op)
            seqpos += len
            refpos += len
        elseif isinsertop(op)
            seqpos += len
        elseif isdeleteop(op)
            refpos += len
        else
            error("operation $(op) is not supported")
        end
        push!(anchors, AlignmentAnchor(seqpos, refpos, op))
    end
    return Alignment(anchors)
end

function Bio.Seq.sequence(rec::BAMRecord)
    seqlen = seqlength(rec)
    data = Vector{UInt64}(cld(seqlen, 16))
    src::Ptr{UInt64} = pointer(rec.data, seqname_length(rec) + n_cigar_op(rec) * 4 + 1)
    for i in 1:endof(data)
        # copy data flipping high and low nybble
        x = unsafe_load(src, i)
        data[i] = (x & 0x0f0f0f0f0f0f0f0f) << 4 | (x & 0xf0f0f0f0f0f0f0f0) >> 4
    end
    return Bio.Seq.DNASequence(data, 1:seqlen, false)
end


"""
    seqlength(rec::BAMRecord)

Return the length of the DNA sequence.
"""
function seqlength(rec)
    return rec.l_seq
end

"""
    qualities(rec::BAMRecord)

Return base qualities of the alignment `rec`.
"""
function qualities(rec::BAMRecord)
    seqlen = seqlength(rec)
    offset = seqname_length(rec) + n_cigar_op(rec) * 4 + cld(seqlen, 2)
    return [reinterpret(Int8, rec.data[i+offset]) for i in 1:seqlen]
end

function Base.getindex(rec::BAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return getvalue(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2]))
end

function Base.setindex!(rec::BAMRecord, val, tag::AbstractString)
    checkkeytag(tag)
    setvalue!(rec.data, auxdata_position(rec), val, UInt8(tag[1]), UInt8(tag[2]))
    return rec
end

function Base.delete!(rec::BAMRecord, tag::AbstractString)
    checkkeytag(tag)
    deletevalue!(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2]))
    return rec
end

function Base.haskey(rec::BAMRecord, tag::AbstractString)
    checkkeytag(tag)
    return findtag(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2])) > 0
end

function optional_fields(rec::BAMRecord)
    return AuxDataDict(rec.data[auxdata_position(rec):data_size(rec)])
end

function auxdata_position(rec)
    seqlen = seqlength(rec)
    return seqname_length(rec) + n_cigar_op(rec) * 4 + cld(seqlen, 2) + seqlen + 1
end

# Return the size of the `.data` field.
function data_size(rec::BAMRecord)
    return rec.block_size - BAM_FIXED_FIELDS_BYTES + sizeof(rec.block_size)
end

# Return the length of alignment.
function alignment_length(rec::BAMRecord)
    offset = seqname_length(rec)
    length::Int = 0
    for i in offset+1:4:offset+n_cigar_op(rec)*4
        x = unsafe_load(Ptr{UInt32}(pointer(rec.data, i)))
        op = Operation(x & 0x0f)
        if ismatchop(op) || isdeleteop(op)
            length += x >> 4
        end
    end
    return length
end

# Return the length of the read name.
function seqname_length(rec)
    return rec.bin_mq_nl & 0xff
end

# Return the number of CIGAR operations.
function n_cigar_op(rec)
    return rec.flag_nc & 0xffff
end
