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

function SAMRecord(data::Vector{UInt8}=UInt8[])
    record = SAMRecord(
        false, data,
        1:0, 1:0, 1:0, 1:0, 1:0,
        1:0, 1:0, 1:0, 1:0, 1:0,
        1:0, UnitRange{Int}[])
    if !isempty(data)
        index!(record)
    end
    return record
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
        println(io, "         sequence: ", sequence(String, record))
        println(io, "   base qualities: ", qualities(String, record))
          print(io, "  optional fields:")
        for field in record.fields
            print(io, ' ', String(record.data[field]))
        end
    else
        print(io, " <not filled>")
    end
end

function Base.print(io::IO, record::SAMRecord)
    write(io, record)
    return nothing
end

function Base.write(io::IO, record::SAMRecord)
    checkfilled(record)
    return unsafe_write(io, pointer(record.data), dataend(record))
end

function isfilled(record::SAMRecord)
    return record.filled
end

function checkfilled(record::SAMRecord)
    if !isfilled(record)
        throw(ArgumentError("unfilled SAM record"))
    end
end

function dataend(record::SAMRecord)
    if record.filled
        if isempty(record.fields)
            return last(record.qual)
        else
            return last(record.fields[end])
        end
    else
        return 0
    end
end

function Base.:(==)(record1::SAMRecord, record2::SAMRecord)
    if record1.filled == record2.filled == false
        return true
    elseif record1.filled == record2.filled == true
        return dataend(record1) == dataend(record2) && memcmp(pointer(record1.data), pointer(record2.data), dataend(record1)) == 0
    else
        return false
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

function rightposition(record::SAMRecord)
    return leftposition(record) + alignment_length(record) - 1
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
    len = length(record.seq)
    if len == 1 && record.data[first(record.seq)] == UInt8('*')
        return Bio.Seq.DNASequence(0)
    else
        ret = Bio.Seq.DNASequence(length(record.seq))
        Bio.Seq.encode_copy!(ret, 1, record.data, first(record.seq), len)
        return ret
    end
end

function sequence(::Type{String}, record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.seq])
end

function qualities(record::SAMRecord)
    checkfilled(record)
    qual = record.data[record.qual]
    for i in 1:endof(qual)
        @inbounds qual[i] -= 33
    end
    return qual
end

function qualities(::Type{String}, record::SAMRecord)
    checkfilled(record)
    return String(record.data[record.qual])
end

function seqlength(record::SAMRecord)
    checkfilled(record)
    if length(record.seq) == 1 && record.data[first(record.seq)] == UInt8('*')
        throw(ArgumentError("no sequence available"))
    end
    return length(record.seq)
end

function alignment(rec::SAMRecord)
    if ismapped(rec)
        return Alignment(cigar(rec), 1, leftposition(rec))
    else
        return Alignment(AlignmentAnchor[])
    end
end

function Base.haskey(record::SAMRecord, tag::AbstractString)
    return findtag(record, tag) > 0
end

function Base.keys(record::SAMRecord)
    checkfilled(record)
    return [String(record.data[first(f):first(f)+1]) for f in record.fields]
end

function Base.values(record::SAMRecord)
    return [record[k] for k in keys(record)]
end

function Base.getindex(record::SAMRecord, tag::AbstractString)
    i = findtag(record, tag)
    if i == 0
        throw(KeyError(tag))
    end
    field = record.fields[i]
    # data type
    typ = record.data[first(field)+3]
    lo = first(field) + 5
    if i == endof(record.fields)
        hi = last(field)
    else
        hi = first(record.fields[i+1]) - 2
    end
    if typ == UInt8('A')
        @assert lo == hi
        return Char(record.data[lo])
    elseif typ == UInt8('i')
        return unsafe_parse_decimal(Int, record.data, lo:hi)
    elseif typ == UInt8('f')
        # TODO: Call a C function directly for speed?
        return parse(Float32, SubString(record.data[lo:hi]))
    elseif typ == UInt8('Z')
        return String(record.data[lo:hi])
    elseif typ == UInt8('H')
        return parse_hexarray(record.data, lo:hi)
    elseif typ == UInt8('B')
        return parse_typedarray(record.data, lo:hi)
    else
        throw(ArgumentError("type code '$(Char(typ))' is not defined"))
    end
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

function optional_fields(record::SAMRecord)
    return Dict(k => record[k] for k in keys(record))
end

# Return the length of alignment.
function alignment_length(record::SAMRecord)
    if length(record.cigar) == 1 && record.data[first(record.cigar)] == UInt8('*')
        return 0
    end
    ret = 0
    len = 0  # operation length
    for i in record.cigar
        c = record.data[i]
        if c ∈ UInt8('0'):UInt8('9')
            len = len * 10 + (c - UInt8('0'))
        else
            op = convert(Operation, Char(c))
            if ismatchop(op) || isdeleteop(op)
                ret += len
                len = 0
            end
        end
    end
    return ret
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

function parse_hexarray(data::Vector{UInt8}, range::UnitRange{Int})
    @assert iseven(length(range))
    ret = Vector{UInt8}(length(range) >> 1)
    byte2hex(b) = b ∈ 0x30:0x39 ? (b - 0x30) : b ∈ 0x41:0x46 ? (b - 0x41 + 0x0A) : error("not in [0-9A-F]")
    j = 1
    for i in first(range):2:last(range)-1
        ret[j] = (byte2hex(data[range[i]]) << 4) | byte2hex(data[range[i+1]])
        j += 1
    end
    return ret
end

function parse_typedarray(data::Vector{UInt8}, range::UnitRange{Int})
    # format: [cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+
    t = data[first(range)]
    xs = split(String(data[first(range)+2:last(range)]))
    if t == UInt8('c')
        return [parse(Int8, x) for x in xs]
    elseif t == UInt8('C')
        return [parse(UInt8, x) for x in xs]
    elseif t == UInt8('s')
        return [parse(Int16, x) for x in xs]
    elseif t == UInt8('S')
        return [parse(UInt16, x) for x in xs]
    elseif t == UInt8('i')
        return [parse(Int32, x) for x in xs]
    elseif t == UInt8('I')
        return [parse(UInt32, x) for x in xs]
    elseif t == UInt8('f')
        return [parse(Float32, x) for x in xs]
    else
        throw(ArgumentError("type code '$(Char(t))' is not defined"))
    end
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end
