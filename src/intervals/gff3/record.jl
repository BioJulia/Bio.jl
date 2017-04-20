# GFF3 Record
# ===========

type Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # :undefiend | :feature | :directive | :comment
    kind::Symbol
    # indexes
    seqid::UnitRange{Int}
    source::UnitRange{Int}
    type_::UnitRange{Int}
    start::UnitRange{Int}
    end_::UnitRange{Int}
    score::UnitRange{Int}
    strand::Int
    phase::Int
    attribute_keys::Vector{UnitRange{Int}}
end

"""
    GFF3.Record()

Create an unfilled GFF3 record.
"""
function Record()
    return Record(
        UInt8[], 1:0, :undefiend,
        # seqid-end_
        1:0, 1:0, 1:0, 1:0, 1:0,
        # score-attribute_keys
        1:0, 0, 0, UnitRange{Int}[])
end

"""
    GFF3.Record(data::Vector{UInt8})

Create a GFF3 record object from `data`.
This function verifies and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(
        data, 1:0, :undefiend,
        # seqid-end_
        1:0, 1:0, 1:0, 1:0, 1:0,
        # score-attribute_keys
        1:0, 0, 0, UnitRange{Int}[])
    index!(record)
    return record
end

"""
    GFF3.Record(str::AbstractString)

Create a GFF3 record object from `str`.
This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return Record(convert(Vector{UInt8}, str))
end

function initialize!(record::Record)
    record.filled = 1:0
    record.kind = :undefiend
    record.seqid = 1:0
    record.source = 1:0
    record.type_ = 1:0
    record.start = 1:0
    record.end_ = 1:0
    record.score = 1:0
    record.strand = 0
    record.phase = 0
    empty!(record.attribute_keys)
    return record
end

function Base.convert(::Type{String}, record::Record)
    return String(record.data[datarange(record)])
end

function Base.convert(::Type{Interval}, record::Record)
    name = Bio.seqname(record)
    lpos = Bio.leftposition(record)
    rpos = Bio.rightposition(record)
    strd = hasstrand(record) ? Bio.Intervals.strand(record) : Bio.Intervals.STRAND_BOTH
    return Interval(name, lpos, rpos, strd, record)
end

function Base.convert(::Type{Interval{Record}}, record::Record)
    return convert(Interval, record)
end

function isfilled(record::Record)
    return !isempty(record.filled)
end

function datarange(record::Record)
    return record.filled
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = datarange(record1)
        r2 = datarange(record2)
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.kind,
        record.seqid,
        record.source,
        record.type_,
        record.start,
        record.end_,
        record.score,
        record.strand,
        record.phase,
        copy(record.attribute_keys))
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        if isfeature(record)
            println(io)
            println(io, "       seqid: ", hasseqid(record) ? seqid(record) : "<missing>")
            println(io, "      source: ", hassource(record) ? source(record) : "<missing>")
            println(io, "        type: ", hasfeaturetype(record) ? featuretype(record) : "<missing>")
            println(io, "       start: ", hasseqstart(record) ? seqstart(record) : "<missing>")
            println(io, "         end: ", hasseqend(record) ? seqend(record) : "<missing>")
            println(io, "       score: ", hasscore(record) ? score(record) : "<missing>")
            println(io, "      strand: ", hasstrand(record) ? strand(record) : "<missing>")
            println(io, "       phase: ", hasphase(record) ? phase(record) : "<missing>")
              print(io, "  attributes:")
            for (key, vals) in attributes(record)
                print(io, ' ', key, '=', join(vals, ','))
            end
        elseif isdirective(record)
            println(io)
              print(io, "   directive: ", String(record))
        elseif iscomment(record)
            println(io)
              print(io, "     comment: ", String(record))
        else
            @assert false
        end
    else
        @assert record.kind == :undefiend
        print(io, " <not filled>")
    end
end


# Accessor functions
# ------------------

"""
    isfeature(record::Record)::Bool

Test if `record` is a feature record.
"""
function isfeature(record::Record)
    return record.kind == :feature
end

"""
    isdirective(record::Record)::Bool

Test if `record` is a directive record.
"""
function isdirective(record::Record)
    return record.kind == :directive
end

"""
    iscomment(record::Record)::Bool

Test if `record` is a comment record.
"""
function iscomment(record::Record)
    return record.kind == :comment
end

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled GFF3 record"))
    end
end

function checkkind(record::Record, kind::Symbol)
    if record.kind != kind
        throw(ArgumentError("$(kind) is expected but $(record.kind) is given"))
    end
end

"""
    seqid(record::Record)::String

Get the sequence ID of `record`.
"""
function seqid(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.seqid)
        missingerror(:seqname)
    end
    return decode(String(record.data[record.seqid]))
end

function hasseqid(record::Record)
    return record.kind == :feature && !ismissing(record, record.seqid)
end

function Bio.seqname(record::Record)
    return seqid(record)
end

function Bio.hasseqname(record::Record)
    return hasseqid(record)
end

"""
    source(record::Record)::String

Get the source of `record`.
"""
function source(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.source)
        missingerror(:source)
    end
    return decode(String(record.data[record.source]))
end

function hassource(record::Record)
    return record.kind == :feature && !ismissing(record, record.source)
end

"""
    featuretype(record::Record)::String

Get the type of `record`.
"""
function featuretype(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.type_)
        missingerror(:featuretype)
    end
    return decode(String(record.data[record.type_]))
end

function hasfeaturetype(record::Record)
    return record.kind == :feature && !ismissing(record, record.type_)
end

"""
    seqstart(record::Record)::Int

Get the start coordinate of `record`.
"""
function seqstart(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.start)
        missingerror(:start)
    end
    return unsafe_parse_decimal(Int, record.data, record.start)
end

function hasseqstart(record::Record)
    return record.kind == :feature && !ismissing(record, record.start)
end

function Bio.leftposition(record::Record)
    return seqstart(record)
end

function Bio.hasleftposition(record::Record)
    return hasseqstart(record)
end

"""
    seqend(record::Record)::Int

Get the end coordinate of `record`.
"""
function seqend(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.end_)
        missingerror(:end)
    end
    return unsafe_parse_decimal(Int, record.data, record.end_)
end

function hasseqend(record::Record)
    return record.kind == :feature && !ismissing(record, record.end_)
end

function Bio.rightposition(record::Record)
    return seqend(record)
end

function Bio.hasrightposition(record::Record)
    return hasseqend(record)
end

"""
    score(record::Record)::Float64

Get the score of `record`
"""
function score(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.score)
        missingerror(:score)
    end
    # TODO: non-copy parsing
    return parse(Float64, String(record.data[record.score]))
end

function hasscore(record::Record)
    return record.kind == :feature && !ismissing(record, record.score)
end

"""
    strand(record::Record)::Bio.Intervals.Strand

Get the strand of `record`.
"""
function strand(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.strand)
        missingerror(:strand)
    end
    return convert(Bio.Intervals.Strand, Char(record.data[record.strand]))
end

function hasstrand(record::Record)
    return record.kind == :feature && !ismissing(record, record.strand)
end

function Bio.Intervals.strand(record::Record)
    return strand(record)
end

"""
    phase(record::Record)::Int

Get the phase of `record`.
"""
function phase(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    if ismissing(record, record.phase)
        missingerror(:phase)
    end
    return unsafe_parse_decimal(Int, record.data, record.phase:record.phase)
end

function hasphase(record::Record)
    return record.kind == :feature && !ismissing(record, record.phase)
end

"""
    attributes(record::Record)::Vector{Pair{String,Vector{String}}}

Get the attributes of `record`.
"""
function attributes(record::Record)
    checkfilled(record)
    checkkind(record, :feature)
    ret = Pair{String,Vector{String}}[]
    for (i, key) in enumerate(record.attribute_keys)
        push!(ret, decode(String(record.data[key])) => attrvals(record, i))
    end
    return ret
end

function attrvals(record::Record, i::Int)
    key = record.attribute_keys[i]
    @assert record.data[last(key)+1] == UInt8('=')
    if i == endof(record.attribute_keys)
        valsrange = last(key)+2:last(datarange(record))
    else
        nextkey = record.attribute_keys[i+1]
        @assert record.data[first(nextkey)-1] == UInt8(';')
        valsrange = last(key)+2:first(nextkey)-2
    end
    stop = last(valsrange)
    ps = first(valsrange)
    pe = searchend(record.data, UInt8(','), ps, stop)
    vals = [decode(String(record.data[ps:pe-1]))]
    while pe ≤ stop
        ps = pe + 1
        pe = searchend(record.data, UInt8(','), ps, stop)
        push!(vals, decode(String(record.data[ps:pe-1])))
    end
    return vals
end

"""
    attributes(record::Record, key::String)::Vector{String}

Get the attributes of `record` with `key`.
"""
function attributes(record::Record, key::String)
    checkfilled(record)
    checkkind(record, :feature)
    for (i, k) in enumerate(record.attribute_keys)
        if isequaldata(key, record.data, k)
            return attrvals(record, i)
        end
    end
    throw(KeyError(key))
end

function searchend(data::Vector{UInt8}, b::UInt8, start::Int, stop::Int)
    p = start
    while p ≤ stop
        if data[p] == b
            return p
        end
        p += 1
    end
    return p
end

"""
    content(record::Record)::String

Get the content of `record`. Leading '#' letters are removed.
"""
function content(record::Record)
    checkfilled(record)
    lo = first(datarange(record))
    hi = last(datarange(record))
    if isfeature(record)
        return decode(String(record.data[lo:hi]))
    elseif isdirective(record)
        return decode(String(record.data[lo+2:hi]))
    elseif iscomment(record)
        return decode(String(record.data[lo+1:hi]))
    else
        @assert false
    end
end

function is_fasta_directive(record::Record)
    return isdirective(record) && isequaldata("##FASTA", record.data, datarange(record))
end

function ismissing(record::Record, index::Int)
    return record.data[index] == UInt8('.')
end

function ismissing(record::Record, range::UnitRange{Int})
    return length(range) == 1 && record.data[first(range)] == UInt8('.')
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

# Decode the percent-encoded string if needed.
function decode(str::String)
    if search(str, '%') > 0
        return URIParser.unescape(str)
    else
        return str
    end
end

# Check if `str == data[range]`
function isequaldata(str::String, data::Vector{UInt8}, range::UnitRange{Int})
    rlen = length(range)
    return rlen == sizeof(str) && memcmp(pointer(data, first(range)), pointer(str), rlen) == 0
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end
