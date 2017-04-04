# BAM Record
# ==========

type Record
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
const FIXED_FIELDS_BYTES = 36

"""
    BAM.Record()

Create an unfilled BAM record.
"""
function Record()
    flag_nc = UInt32(SAM.FLAG_UNMAP) << 16
    return Record(0, -1, -1, 0, flag_nc, 0, -1, -1, 0, UInt8[], String[])
end

function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record()
    unsafe_copy!(pointer_from_objref(record), pointer(data), FIXED_FIELDS_BYTES)
    dsize = data_size(record)
    resize!(record.data, dsize)
    unsafe_copy!(pointer(record.data), pointer(data, FIXED_FIELDS_BYTES + 1), dsize)
    return record
end

# Return the size of the `.data` field.
function data_size(record::Record)
    return record.block_size - FIXED_FIELDS_BYTES + sizeof(record.block_size)
end

function Bio.isfilled(record::Record)
    return record.block_size != 0
end

# NOTE: this does not copy `refseqnames`.
function Base.copy(rec::Record)
    return Record(
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

function Base.isless(rec1::Record, rec2::Record)
    # compared by left-most position of an alignment
    if rec1.refid == rec2.refid
        return isless(rec1.pos, rec2.pos)
    else
        return isless(rec1.refid, rec2.refid)
    end
end

function Base.:(==)(rec1::Record, rec2::Record)
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

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        # TODO
    else
        print(io, " <not filled>")
    end
end

"""
    ismapped(record::Record)

Test if `record` is mapped.
"""
function ismapped(record::Record)
    return flag(record) & SAM.FLAG_UNMAP == 0
end

"""
    refid(record::Record)::Int

Get the ID of a reference sequence of `record`.

The ID is 1-based (i.e. the first sequence is 1) and is 0 for a record without a mapping position.
"""
function refid(record::Record)::Int
    checkfilled(record)
    return record.refid + 1
end

# TODO: rname?

"""
    refindex(rec::Record)

Return the index of a reference sequence that `rec` is mapped onto.

The index is 1-based and will be 0 for an alignment without mapping position.
"""
function refindex(rec::Record)
    return rec.refid + 1
end

function nextrefindex(rec::Record)
    return rec.next_refid + 1
end

"""
    refname(rec::Record)

Return the name of a reference sequence that `rec` is mapped onto.

If `rec` is unmapped, it returns `"*"` like SAM records.
"""
function refname(rec::Record)
    i = refindex(rec)
    if i == 0
        return "*"
    else
        return rec.refseqnames[i]
    end
end

function nextrefname(rec::Record)
    i = nextrefindex(rec)
    if i == 0
        return "*"
    else
        return rec.refseqnames[i]
    end
end

"""
    pos(record::Record)::Int

Get the 1-based leftmost mapping position of `record`.
"""
function pos(record::Record)::Int
    checkfilled(record)
    return record.pos + 1
end

function haspos(record::Record)
    return isfilled(record)
end

function Bio.leftposition(record::Record)
    return pos(record)
end

function Bio.hasleftposition(record::Record)
    return haspos(record)
end

"""
    rightposition(rec::Record)

Return the rightmost mapping position of `rec`.
"""
function Bio.rightposition(record::Record)
    checkfilled(record)
    return Int32(pos(record) + alignment_length(record) - 1)
end

function nextleftposition(rec::Record)
    return rec.next_pos + 1
end

"""
    mapq(record::Record)::UInt8

Get the mapping quality of `record`.
"""
function mapq(record::Record)::UInt8
    return UInt8((record.bin_mq_nl >> 8) & 0xff)
end

function hasmapq(record::Record)
    return isfilled(record)
end

"""
    flag(record::Record)::UInt16

Get the bitwise flag of `record`.
"""
function flag(record::Record)::UInt16
    checkfilled(record)
    return UInt16(record.flag_nc >> 16)
end

function hasflag(record::Record)
    return isfilled(record)
end

"""
    tlen(record::Record)

Get the template length of `record`.
"""
function tlen(record::Record)
    checkfilled(record)
    return record.tlen
end

function hastlen(record::Record)
    return isfilled(record)
end

# TODO: is this needed?
function Bio.seqname(rec::Record)
    # drop the last NUL character
    return unsafe_string(pointer(rec.data), max(seqname_length(rec) - 1, 0))
end

"""
    cigar_rle(record::Record)::Tuple{Vector{Bio.Align.Operation},Vector{Int}}

Get a run-length encoded tuple `(ops, lens)` of the CIGAR string in `record`.

See also `BAM.cigar`.
"""
function cigar_rle(record::Record)
    checkfilled(record)
    offset = seqname_length(record)
    ops = Bio.Align.Operation[]
    lens = Int[]
    for i in offset+1:4:offset+n_cigar_op(record)*4
        x = unsafe_load(Ptr{UInt32}(pointer(record.data, i)))
        op = Bio.Align.Operation(x & 0x0f)
        push!(ops, op)
        push!(lens, x >> 4)
    end
    return ops, lens
end

"""
    cigar(record::Record)::String

Get the CIGAR string of `record`.

See also `cigar_rle`.
"""
function cigar(record::Record)::String
    buf = IOBuffer()
    for (op, len) in zip(cigar_rle(record)...)
        print(buf, len, Char(op))
    end
    return takebuf_string(buf)
end

"""
    alignment(record::Record)::Bio.Align.Alignment

Get the alignment of `record`.
"""
function alignment(record::Record)::Bio.Align.Alignment
    checkfilled(record)
    if !ismapped(record)
        return Bio.Align.Alignment(Bio.Align.AlignmentAnchor[])
    end
    seqpos = 0
    refpos = pos(record) - 1
    anchors = [Bio.Align.AlignmentAnchor(seqpos, refpos, Bio.Align.OP_START)]
    for (op, len) in zip(cigar_rle(record)...)
        if Bio.Align.ismatchop(op)
            seqpos += len
            refpos += len
        elseif Bio.Align.isinsertop(op)
            seqpos += len
        elseif Bio.Align.isdeleteop(op)
            refpos += len
        else
            error("operation $(op) is not supported")
        end
        push!(anchors, Bio.Align.AlignmentAnchor(seqpos, refpos, op))
    end
    return Bio.Align.Alignment(anchors)
end

"""
    seq(record::Record)::Bio.Seq.DNASequence

Get the segment sequence of `record`.
"""
function seq(record::Record)
    checkfilled(record)
    seqlen = seqlength(record)
    data = Vector{UInt64}(cld(seqlen, 16))
    src::Ptr{UInt64} = pointer(record.data, seqname_length(record) + n_cigar_op(record) * 4 + 1)
    for i in 1:endof(data)
        # copy data flipping high and low nybble
        x = unsafe_load(src, i)
        data[i] = (x & 0x0f0f0f0f0f0f0f0f) << 4 | (x & 0xf0f0f0f0f0f0f0f0) >> 4
    end
    return Bio.Seq.DNASequence(data, 1:seqlen, false)
end

function hasseq(record::Record)
    return isfilled(record)
end

function Bio.sequence(record::Record)
    return seq(record)
end

function Bio.hassequence(record::Record)
    return hasseq(record)
end

"""
    seqlength(record::Record)::Int

Get the sequence length of `record`.
"""
function seqlength(record::Record)::Int
    checkfilled(record)
    return record.l_seq % Int
end

"""
    qual(record::Record)::Vector{UInt8}

Get the base quality of  `record`.
"""
function qual(record::Record)::Vector{UInt8}
    seqlen = seqlength(record)
    offset = seqname_length(record) + n_cigar_op(record) * 4 + cld(seqlen, 2)
    return [reinterpret(Int8, record.data[i+offset]) for i in 1:seqlen]
end

function hasqual(record::Record)
    return isfilled(record)
end

function Base.getindex(rec::Record, tag::AbstractString)
    checkkeytag(tag)
    return getvalue(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2]))
end

function Base.setindex!(rec::Record, val, tag::AbstractString)
    checkkeytag(tag)
    setvalue!(rec.data, auxdata_position(rec), val, UInt8(tag[1]), UInt8(tag[2]))
    return rec
end

function Base.delete!(rec::Record, tag::AbstractString)
    checkkeytag(tag)
    deletevalue!(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2]))
    return rec
end

function Base.haskey(rec::Record, tag::AbstractString)
    checkkeytag(tag)
    return findtag(rec.data, auxdata_position(rec), UInt8(tag[1]), UInt8(tag[2])) > 0
end

function optional_fields(rec::Record)
    return AuxDataDict(rec.data[auxdata_position(rec):data_size(rec)])
end

function auxdata_position(rec)
    seqlen = seqlength(rec)
    return seqname_length(rec) + n_cigar_op(rec) * 4 + cld(seqlen, 2) + seqlen + 1
end

# Return the length of alignment.
function alignment_length(rec::Record)
    offset = seqname_length(rec)
    length::Int = 0
    for i in offset+1:4:offset+n_cigar_op(rec)*4
        x = unsafe_load(Ptr{UInt32}(pointer(rec.data, i)))
        op = Bio.Align.Operation(x & 0x0f)
        if Bio.Align.ismatchop(op) || Bio.Align.isdeleteop(op)
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

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled BAM record"))
    end
end
