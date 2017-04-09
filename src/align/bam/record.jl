# BAM Record
# ==========

"""
    BAM.Record()

Create an unfilled BAM record.
"""
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
    reader::Reader

    function Record()
        return new(0, 0, 0, 0, 0, 0, 0, 0, 0, UInt8[])
    end
end

# the data size of fixed-length fields (block_size-tlen)
const FIXED_FIELDS_BYTES = 36

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

function Base.copy(record::Record)
    copy = Record()
    copy.block_size = record.block_size
    copy.refid      = record.refid
    copy.pos        = record.pos
    copy.bin_mq_nl  = record.bin_mq_nl
    copy.flag_nc    = record.flag_nc
    copy.l_seq      = record.l_seq
    copy.next_refid = record.next_refid
    copy.next_pos   = record.next_pos
    copy.tlen       = record.tlen
    copy.data       = record.data[1:data_size(record)]
    if isdefined(record, :reader)
        copy.reader = record.reader
    end
    return copy
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "      template name: ", tempname(record))
        println(io, "               flag: ", flag(record))
        println(io, "       reference ID: ", refid(record))
        println(io, "           position: ", position(record))
        println(io, "    mapping quality: ", mappingquality(record))
        println(io, "              CIGAR: ", cigar(record))
        println(io, "  next reference ID: ", nextrefid(record))
        println(io, "      next position: ", nextposition(record))
        println(io, "    template length: ", templength(record))
        println(io, "           sequence: ", sequence(record))
        # TODO: pretty print base quality
        println(io, "       base quality: ", quality(record))
          print(io, "     auxiliary data:")
        for field in keys(auxdata(record))
            print(io, ' ', field, '=', record[field])
        end
    else
        print(io, " <not filled>")
    end
end

function Base.read!(reader::Reader, record::Record)
    return _read!(reader, record)
end


# Accessor Fuctions
# -----------------

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
    ismapped(record::Record)::Bool

Test if `record` is mapped.
"""
function ismapped(record::Record)
    return flag(record) & SAM.FLAG_UNMAP == 0
end

"""
    refid(record::Record)::Int

Get the reference sequence ID of `record`.

The ID is 1-based (i.e. the first sequence is 1) and is 0 for a record without a mapping position.

See also: `BAM.rname`
"""
function refid(record::Record)::Int
    checkfilled(record)
    return record.refid + 1
end

function hasrefid(record::Record)
    return isfilled(record)
end

"""
    refname(record::Record)::String

Get the reference sequence name of `record`.

See also: `BAM.refid`
"""
function refname(record::Record)::String
    checkfilled(record)
    id = refid(record)
    if id == 0
        throw(ArgumentError("record is not mapped"))
    elseif !isdefined(record, :reader)
        throw(ArgumentError("reader is not defined"))
    end
    return record.reader.refseqnames[id]
end

function hasrefname(record::Record)
    return hasrefid(record)
end

"""
    position(record::Record)::Int

Get the 1-based leftmost mapping position of `record`.
"""
function position(record::Record)::Int
    checkfilled(record)
    return record.pos + 1
end

function hasposition(record::Record)
    return isfilled(record)
end

"""
    rightposition(record::Record)::Int

Get the 1-based rightmost mapping position of `record`.
"""
function rightposition(record::Record)::Int
    checkfilled(record)
    return Int32(position(record) + alignlength(record) - 1)
end

function hasrightposition(record::Record)
    return isfilled(record) && ismapped(record)
end

"""
    isnextmapped(record::Record)::Bool

Test if the mate/next read of `record` is mapped.
"""
function isnextmapped(record::Record)::Bool
    return isfilled(record) && (flag(record) & FLAG_MUNMAP == 0)
end

"""
    nextrefid(record::Record)::Int

Get the next/mate reference sequence ID of `record`.
"""
function nextrefid(record::Record)::Int
    checkfilled(record)
    return record.next_refid + 1
end

function hasnextrefid(record::Record)
    return isfilled(record)
end

"""
    nextrefname(record::Record)::String

Get the reference name of the mate/next read of `record`.
"""
function nextrefname(record::Record)::String
    checkfilled(record)
    id = nextrefid(record)
    if id == 0
        throw(ArgumentError("next record is not mapped"))
    elseif !isdefined(record, :reader)
        throw(ArgumentError("reader is not defined"))
    end
    return record.reader.refseqnames[id]
end

function hasnextrefname(record::Record)
    return isfilled(record) && isnextmapped(record)
end

"""
    nextposition(record::Record)::Int

Get the 1-based leftmost mapping position of the next/mate read of `record`.
"""
function nextposition(record::Record)::Int
    checkfilled(record)
    return record.next_pos + 1
end

function hasnextposition(record::Record)
    return isfilled(record)
end

"""
    mappingquality(record::Record)::UInt8

Get the mapping quality of `record`.
"""
function mappingquality(record::Record)::UInt8
    return UInt8((record.bin_mq_nl >> 8) & 0xff)
end

function hasmappingquality(record::Record)
    return isfilled(record)
end

"""
    cigar(record::Record)::String

Get the CIGAR string of `record`.

See also `BAM.cigar_rle`.
"""
function cigar(record::Record)::String
    buf = IOBuffer()
    for (op, len) in zip(cigar_rle(record)...)
        print(buf, len, Char(op))
    end
    return takebuf_string(buf)
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
    alignment(record::Record)::Bio.Align.Alignment

Get the alignment of `record`.
"""
function alignment(record::Record)::Bio.Align.Alignment
    checkfilled(record)
    if !ismapped(record)
        return Bio.Align.Alignment(Bio.Align.AlignmentAnchor[])
    end
    seqpos = 0
    refpos = position(record) - 1
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

function hasalignment(record::Record)
    return ismapped(record)
end

"""
    alignlength(record::Record)::Int

Get the alignment length of `record`.
"""
function alignlength(record::Record)::Int
    offset = seqname_length(record)
    length::Int = 0
    for i in offset+1:4:offset+n_cigar_op(record)*4
        x = unsafe_load(Ptr{UInt32}(pointer(record.data, i)))
        op = Bio.Align.Operation(x & 0x0f)
        if Bio.Align.ismatchop(op) || Bio.Align.isdeleteop(op)
            length += x >> 4
        end
    end
    return length
end

"""
    tempname(record::Record)::String

Get the query template name of `record`.
"""
function tempname(record::Record)::String
    checkfilled(record)
    # drop the last NUL character
    return unsafe_string(pointer(record.data), max(seqname_length(record) - 1, 0))
end

function hastempname(record::Record)
    return isfilled(record)
end

"""
    templength(record::Record)::Int

Get the template length of `record`.
"""
function templength(record::Record)::Int
    checkfilled(record)
    return record.tlen
end

function hastemplength(record::Record)
    return isfilled(record)
end

"""
    sequence(record::Record)::Bio.Seq.DNASequence

Get the segment sequence of `record`.
"""
function sequence(record::Record)::Bio.Seq.DNASequence
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

function hassequence(record::Record)
    return isfilled(record)
end

"""
    seqlength(record::Record)::Int

Get the sequence length of `record`.
"""
function seqlength(record::Record)::Int
    checkfilled(record)
    return record.l_seq % Int
end

function hasseqlength(record::Record)
    return isfilled(record)
end

"""
    quality(record::Record)::Vector{UInt8}

Get the base quality of  `record`.
"""
function quality(record::Record)::Vector{UInt8}
    checkfilled(record)
    seqlen = seqlength(record)
    offset = seqname_length(record) + n_cigar_op(record) * 4 + cld(seqlen, 2)
    return [reinterpret(Int8, record.data[i+offset]) for i in 1:seqlen]
end

function hasquality(record::Record)
    return isfilled(record)
end

"""
    auxdata(record::Record)::BAM.AuxData

Get the auxiliary data of `record`.
"""
function auxdata(record::Record)
    checkfilled(record)
    return AuxData(record.data[auxdata_position(record):data_size(record)])
end

function hasauxdata(record::Record)
    return isfilled(record)
end

function Base.getindex(record::Record, tag::AbstractString)
    checkauxtag(tag)
    return getauxvalue(record.data, auxdata_position(record), UInt8(tag[1]), UInt8(tag[2]))
end

function Base.haskey(record::Record, tag::AbstractString)
    checkauxtag(tag)
    return findauxtag(record.data, auxdata_position(record), UInt8(tag[1]), UInt8(tag[2])) > 0
end

function Base.keys(record::Record)
    return collect(keys(auxdata(record)))
end

function Base.values(record::Record)
    return [record[key] for key in keys(record)]
end


# Bio Methods
# -----------

function Bio.isfilled(record::Record)
    return record.block_size != 0
end

function Bio.seqname(record::Record)
    return tempname(record)
end

function Bio.hasseqname(record::Record)
    return hastempname(record)
end

function Bio.sequence(record::Record)
    return sequence(record)
end

function Bio.hassequence(record::Record)
    return hassequence(record)
end

function Bio.leftposition(record::Record)
    return position(record)
end

function Bio.hasleftposition(record::Record)
    return hasposition(record)
end

function Bio.rightposition(record::Record)
    return rightposition(record)
end

function Bio.hasrightposition(record::Record)
    return hasrightposition(record)
end


# Helper Functions
# ----------------

# Return the size of the `.data` field.
function data_size(record::Record)
    if isfilled(record)
        return record.block_size - FIXED_FIELDS_BYTES + sizeof(record.block_size)
    else
        return 0
    end
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled BAM record"))
    end
end

function auxdata_position(record::Record)
    seqlen = seqlength(record)
    return seqname_length(record) + n_cigar_op(record) * 4 + cld(seqlen, 2) + seqlen + 1
end

# Return the length of the read name.
function seqname_length(record::Record)
    return record.bin_mq_nl & 0xff
end

# Return the number of CIGAR operations.
function n_cigar_op(record::Record)
    return record.flag_nc & 0xffff
end
