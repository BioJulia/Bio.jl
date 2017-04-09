# SAM Record
# ==========

type Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
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

"""
    SAM.Record()

Create an unfilled SAM record.
"""
function Record()
    return Record(
        UInt8[], 1:0,
        # qname-mapq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # cigar-seq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # qual and fields
        1:0, UnitRange{Int}[])
end

"""
    SAM.Record(data::Vector{UInt8})

Create a SAM record from `data`.
This function verifies the format and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(
        data, 1:0,
        # qname-mapq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # cigar-seq
        1:0, 1:0, 1:0, 1:0, 1:0,
        # qual and fields
        1:0, UnitRange{Int}[])
    index!(record)
    return record
end

"""
    SAM.Record(str::AbstractString)

Create a SAM record from `str`.
This function verifies the format and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return Record(convert(Vector{UInt8}, str))
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "    template name: ", hastempname(record) ? tempname(record) : "<missing>")
        println(io, "             flag: ", hasflag(record) ? flag(record) : "<missing>")
        println(io, "        reference: ", hasrefname(record) ? refname(record) : "<missing>")
        println(io, "         position: ", hasposition(record) ? position(record) : "<missing>")
        println(io, "  mapping quality: ", hasmappingquality(record) ? mappingquality(record) : "<missing>")
        println(io, "            CIGAR: ", hascigar(record) ? cigar(record) : "<missing>")
        println(io, "   next reference: ", hasnextrefname(record) ? nextrefname(record) : "<missing>")
        println(io, "    next position: ", hasnextposition(record) ? nextposition(record) : "<missing>")
        println(io, "  template length: ", hastemplength(record) ? templength(record) : "<missing>")
        println(io, "         sequence: ", hassequence(record) ? sequence(String, record) : "<missing>")
        println(io, "     base quality: ", hasquality(record) ? quality(String, record) : "<missing>")
          print(io, "   auxiliary data:")
        for field in record.fields
            print(io, ' ', String(record.data[field]))
        end
    else
        print(io, " <not filled>")
    end
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.write(io::IO, record::Record)
    checkfilled(record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.copy(record::Record)
    return Record(
        copy(record.data),
        record.filled,
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


# Accessor Functions
# ------------------

"""
    flag(record::Record)::UInt16

Get the bitwise flag of `record`.
"""
function flag(record::Record)::UInt16
    checkfilled(record)
    return unsafe_parse_decimal(UInt16, record.data, record.flag)
end

function hasflag(record::Record)
    return isfilled(record)
end

"""
    ismapped(record::Record)::Bool

Test if `record` is mapped.
"""
function ismapped(record::Record)::Bool
    return isfilled(record) && (flag(record) & FLAG_UNMAP == 0)
end

"""
    refname(record::Record)::String

Get the reference sequence name of `record`.
"""
function refname(record::Record)
    checkfilled(record)
    if ismissing(record, record.rname)
        missingerror(:refname)
    end
    return String(record.data[record.rname])
end

function hasrefname(record::Record)
    return isfilled(record) && !ismissing(record, record.rname)
end

"""
    position(record::Record)::Int

Get the 1-based leftmost mapping position of `record`.
"""
function position(record::Record)::Int
    checkfilled(record)
    pos = unsafe_parse_decimal(Int, record.data, record.pos)
    if pos == 0
        missingerror(:position)
    end
    return pos
end

function hasposition(record::Record)
    return isfilled(record) && (length(record.pos) != 1 || record.data[first(record.pos)] != UInt8('0'))
end

"""
    rightposition(record::Record)::Int

Get the 1-based rightmost mapping position of `record`.
"""
function rightposition(record::Record)
    return position(record) + alignlength(record) - 1
end

function hasrightposition(record::Record)
    return hasposition(record) && hasalignment(record)
end

"""
    isnextmapped(record::Record)::Bool

Test if the mate/next read of `record` is mapped.
"""
function isnextmapped(record::Record)::Bool
    return isfilled(record) && (flag(record) & FLAG_MUNMAP == 0)
end

"""
    nextrefname(record::Record)::String

Get the reference name of the mate/next read of `record`.
"""
function nextrefname(record::Record)::String
    checkfilled(record)
    if ismissing(record, record.rnext)
        missingerror(:nextrefname)
    end
    return String(record.data[record.rnext])
end

function hasnextrefname(record::Record)
    return isfilled(record) && !ismissing(record, record.rnext)
end

"""
    nextposition(record::Record)::Int

Get the position of the mate/next read of `record`.
"""
function nextposition(record::Record)::Int
    checkfilled(record)
    pos = unsafe_parse_decimal(Int, record.data, record.pnext)
    if pos == 0
        missingerror(:nextposition)
    end
    return pos
end

function hasnextposition(record::Record)
    return isfilled(record) && (length(record.pnext) != 1 || record.data[first(record.pnext)] != UInt8('0'))
end

"""
    mappingquality(record::Record)::UInt8

Get the mapping quality of `record`.
"""
function mappingquality(record::Record)::UInt8
    checkfilled(record)
    qual = unsafe_parse_decimal(UInt8, record.data, record.mapq)
    if qual == 0xff
        missingerror(:mappingquality)
    end
    return qual
end

function hasmappingquality(record::Record)
    return isfilled(record) && unsafe_parse_decimal(UInt8, record.data, record.mapq) != 0xff
end

"""
    cigar(record::Record)::String

Get the CIGAR string of `record`.
"""
function cigar(record::Record)::String
    checkfilled(record)
    if ismissing(record, record.cigar)
        missingerror(:cigar)
    end
    return String(record.data[record.cigar])
end

function hascigar(record::Record)
    return isfilled(record) && !ismissing(record, record.cigar)
end

"""
    alignment(record::Record)::Bio.Align.Alignment

Get the alignment of `record`.
"""
function alignment(record::Record)::Bio.Align.Alignment
    if ismapped(record)
        return Bio.Align.Alignment(cigar(record), 1, position(record))
    else
        return Bio.Align.Alignment(Bio.Align.AlignmentAnchor[])
    end
end

function hasalignment(record::Record)
    return isfilled(record) && hascigar(record)
end

"""
    alignlength(record::Record)::Int

Get the alignment length of `record`.
"""
function alignlength(record::Record)::Int
    if length(record.cigar) == 1 && record.data[first(record.cigar)] == UInt8('*')
        return 0
    end
    ret::Int = 0
    len = 0  # operation length
    for i in record.cigar
        c = record.data[i]
        if c ∈ UInt8('0'):UInt8('9')
            len = len * 10 + (c - UInt8('0'))
        else
            op = convert(Bio.Align.Operation, Char(c))
            if Bio.Align.ismatchop(op) || Bio.Align.isdeleteop(op)
                ret += len
                len = 0
            end
        end
    end
    return ret
end

"""
    tempname(record::Record)::String

Get the query template name of `record`.
"""
function tempname(record::Record)::String
    checkfilled(record)
    if ismissing(record, record.qname)
        missingerror(:tempname)
    end
    return String(record.data[record.qname])
end

function hastempname(record::Record)
    return isfilled(record) && !ismissing(record, record.qname)
end

"""
    templength(record::Record)::Int

Get the template length of `record`.
"""
function templength(record::Record)::Int
    checkfilled(record)
    len = unsafe_parse_decimal(Int, record.data, record.tlen)
    if len == 0
        missingerror(:tlen)
    end
    return len
end

function hastemplength(record::Record)
    return isfilled(record) && (length(record.tlen) != 1 || record.data[first(record.tlen)] != UInt8('0'))
end

"""
    sequence(record::Record)::Bio.Seq.DNASequence

Get the segment sequence of `record`.
"""
function sequence(record::Record)::Bio.Seq.DNASequence
    checkfilled(record)
    if ismissing(record, record.seq)
        missingerror(:sequence)
    end
    seqlen = length(record.seq)
    ret = Bio.Seq.DNASequence(seqlen)
    Bio.Seq.encode_copy!(ret, 1, record.data, first(record.seq), seqlen)
    return ret
end

function hassequence(record::Record)
    return isfilled(record) && !ismissing(record, record.seq)
end
"""
    sequence(::Type{String}, record::Record)::String

Get the segment sequence of `record` as `String`.
"""
function sequence(::Type{String}, record::Record)::String
    checkfilled(record)
    return String(record.data[record.seq])
end

"""
    seqlength(record::Record)::Int

Get the sequence length of `record`.
"""
function seqlength(record::Record)::Int
    checkfilled(record)
    if ismissing(record, record.seq)
        missingerror(:seq)
    end
    return length(record.seq)
end

function hasseqlength(record::Record)
    return isfilled(record) && !ismissing(record, record.seq)
end

"""
    quality(record::Record)::Vector{UInt8}

Get the Phred-scaled base quality of `record`.
"""
function quality(record::Record)::Vector{UInt8}
    checkfilled(record)
    if ismissing(record, record.qual)
        missingerror(:quality)
    end
    qual = record.data[record.qual]
    for i in 1:endof(qual)
        @inbounds qual[i] -= 33
    end
    return qual
end

function hasquality(record::Record)
    return isfilled(record) && !ismissing(record, record.qual)
end

"""
    quality(::Type{String}, record::Record)::String

Get the ASCII-encoded base quality of `record`.
"""
function quality(::Type{String}, record::Record)::String
    checkfilled(record)
    return String(record.data[record.qual])
end

"""
    auxdata(record::Record)::Dict{String,Any}

Get the auxiliary data (optional fields) of `record`.
"""
function auxdata(record::Record)::Dict{String,Any}
    checkfilled(record)
    return Dict(k => record[k] for k in keys(record))
end

function Base.haskey(record::Record, tag::AbstractString)
    return findauxtag(record, tag) > 0
end

function Base.getindex(record::Record, tag::AbstractString)
    i = findauxtag(record, tag)
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

function Base.keys(record::Record)
    checkfilled(record)
    return [String(record.data[first(f):first(f)+1]) for f in record.fields]
end

function Base.values(record::Record)
    return [record[k] for k in keys(record)]
end


# Bio Methods
# -----------

function Bio.isfilled(record::Record)
    return !isempty(record.filled)
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

function Bio.rightposition(record::Record)
    return rightposition(record)
end

function Bio.hasrightposition(record::Record)
    return hasrightposition(record)
end

function Bio.leftposition(record::Record)
    return position(record)
end

function Bio.hasleftposition(record::Record)
    return hasposition(record)
end


# Helper Functions
# ----------------

function initialize!(record::Record)
    record.filled = 1:0
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

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled SAM record"))
    end
end

function findauxtag(record::Record, tag::AbstractString)
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

function ismissing(record::Record, range::UnitRange{Int})
    return length(range) == 1 && record.data[first(range)] == UInt8('*')
end
