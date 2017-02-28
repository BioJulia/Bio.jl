# SAM Record
# ==========

type SAMRecord
    filled::Bool
    data::Vector{UInt8}
    qname::UnitRange{Int}
    flag::UnitRange{Int}
    rname::UnitRange{Int}
    pos::UnitRange{Int}
    mapq::UnitRange{Int}
    cigar::UnitRange{Int}
    rnext::UnitRange{Int}
    pnext::UnitRange{Int}
    tlen::UnitRange{Int}
    seq::UnitRange{Int}
    qual::UnitRange{Int}
    fields::Vector{UnitRange{Int}}
end

function SAMRecord()
    return SAMRecord(
        false,
        UInt8[],
        1:0, 1:0, 1:0, 1:0, 1:0,
        1:0, 1:0, 1:0, 1:0, 1:0,
        1:0, UnitRange{Int}[])
end

function initialize!(record::SAMRecord)
    record.filled = false
    record.qname = 1:0
    record.flag = 1:0
    record.rname = 1:0
    record.pos = 1:0
    record.mapq = 1:0
    record.cigar = 1:0
    record.rnext = 1:0
    record.pnext = 1:0
    record.tlen = 1:0
    record.seq = 1:0
    record.qual = 1:0
    empty!(record.fields)
    return record
end

function Base.show(io::IO, record::SAMRecord)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "    sequence name: ", seqname(record))
        println(io, "             flag: ", flag(record))
        println(io, "        reference: ", refname(record))
        println(io, "         position: ", leftposition(record))
        println(io, "  mapping quality: ", mappingquality(record))
        println(io, "            CIGAR: ", cigar(record))
        println(io, "   next reference: ", nextrefname(record))
        println(io, "    next position: ", nextleftposition(record))
        println(io, "  template length: ", templatelength(record))
        println(io, "         sequence: ", sequence(record))
        println(io, "   base qualities: ", qualities(record))
          print(io, "  optional fields:")
        for field in record.fields
            print(io, ' ', String(record.data[field]))
        end
    else
        print(io, " <not filled>")
    end
end

function isfilled(record::SAMRecord)
    return record.filled
end

function checkfilled(record::SAMRecord)
    if !isfilled(record)
        throw(ArgumentError("unfilled SAM record"))
    end
end

# TODO
function Base.isless(rec1::SAMRecord, rec2::SAMRecord)
    # compared by the left-most position of an alignment
    if rec1.name == rec2.name
        return isless(rec1.pos, rec2.pos)
    else
        return isless(rec1.name, rec2.name)
    end
end

function Base.copy(record::SAMRecord)
    return SAMRecord(
        record.filled,
        copy(record.data),
        record.qname,
        record.flag,
        record.rname,
        record.pos,
        record.mapq,
        record.cigar,
        record.rnext,
        record.pnext,
        record.tlen,
        record.seq,
        record.qual,
        copy(record.fields))
end

function ismapped(rec::SAMRecord)
    return flag(rec) & SAM_FLAG_UNMAP == 0
end

function Bio.seqname(record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.qname])
end

function flag(record::SAMRecord)
    checkfilled(record)
    return unsafe_parse_decimal(UInt16, record.data, record.flag)
end

function refname(record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.rname])
end

function leftposition(record::SAMRecord)
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.pos)
end

function mappingquality(record::SAMRecord)
    checkfilled(record)
    return unsafe_parse_decimal(UInt8, record.data, record.mapq)
end

function cigar(record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.cigar])
end

function nextrefname(record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.rnext])
end

function nextleftposition(record::SAMRecord)
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.pnext)
end

function templatelength(record::SAMRecord)
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.tlen)
end

function sequence(record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.seq])
end

function qualities(record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.qual])
end

function seqlength(record::SAMRecord)
    checkfilled(record)
    if length(record.seq) == 1 && record.data[first(record.seq)] == UInt8('*')
        throw(ArgumentError("no sequence available"))
    end
    return length(rec.seq)
end

function alignment(rec::SAMRecord)
    if ismapped(rec)
        return Alignment(cigar(rec), 1, leftposition(rec))
    else
        return Alignment(AlignmentAnchor[])
    end
end

function Base.keys(record::SAMRecord)
    checkfilled(record)
    return [String(record.data[first(f):first(f)+1]) for f in record.fields]
end

# TODO: values
function Base.values(record::SAMRecord)
end

function Base.haskey(record::SAMRecord, tag::AbstractString)
    return findtag(record, tag) > 0
end

function Base.getindex(record::SAMRecord, tag::AbstractString)
    i = findtag(record, tag)
    if i == 0
        throw(KeyError(tag))
    end
    field = record.fields[i]
    # TODO: type
    #typ = record.data[first(field)+3]
    lo = first(field) + 5
    if i == endof(record.fields)
        hi = last(field)
    else
        hi = first(record.fields[i+1]) - 2
    end
    return String(record.data[lo:hi])
end

function findtag(record::SAMRecord, tag::AbstractString)
    checkfilled(record)
    if sizeof(tag) != 2
        return 0
    end
    t1, t2 = UInt8(tag[1]), UInt8(tag[2])
    for (i, field) in enumerate(record.fields)
        p = first(field)
        if record.data[p] == t1 && record.data[p+1] == t2
            return i
        end
    end
    return 0
end

# TODO
function optional_fields(rec::SAMRecord)
    return rec.optional_fields
end

# TODO
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

# r"[0-9]+" must match `data[range]`.
function unsafe_parse_decimal{T<:Unsigned}(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int})
    x = zero(T)
    @inbounds for i in range
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return x
end

# r"[-+]?[0-9]+" must match `data[range]`.
function unsafe_parse_decimal{T<:Signed}(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int})
    lo = first(range)
    if data[lo] == UInt8('-')
        sign = T(-1)
        lo += 1
    elseif data[lo] == UInt8('+')
        sign = T(+1)
        lo += 1
    else
        sign = T(+1)
    end
    x = zero(T)
    @inbounds for i in lo:last(range)
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return sign * x
end

