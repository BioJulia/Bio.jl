# FASTA Record
# ============

type Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
    identifier::UnitRange{Int}
    description::UnitRange{Int}
    sequence::UnitRange{Int}
end

"""
    FASTA.Record()

Create an unfilled FASTA record.
"""
function Record()
    return Record(UInt8[], 1:0, 1:0, 1:0, 1:0)
end

"""
    FASTA.Record(data::Vector{UInt8})

Create a FASTA record object from `data`.

This function verifies and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(data, 1:0, 1:0, 1:0, 1:0)
    index!(record)
    return record
end

"""
    FASTA.Record(str::AbstractString)

Create a FASTA record object from `str`.

This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return Record(convert(Vector{UInt8}, str))
end

function Base.convert(::Type{String}, record::Record)
    return String(record.data[record.filled])
end

"""
    FASTA.Record(identifier, sequence)

Create a FASTA record object from `identifier` and `sequence`.
"""
function Record(identifier::AbstractString, sequence)
    return Record(identifier, nothing, sequence)
end

"""
    FASTA.Record(identifier, description, sequence)

Create a FASTA record object from `identifier`, `description` and `sequence`.
"""
function Record(identifier::AbstractString, description::Union{AbstractString,Void}, sequence)
    buf = IOBuffer()
    print(buf, '>', strip(identifier))
    if description != nothing
        print(buf, ' ', description)
    end
    print(buf, '\n')
    print(buf, sequence, '\n')
    return Record(takebuf_array(buf))
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = record1.filled
        r2 = record2.filled
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.identifier,
        record.description,
        record.sequence)
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled), length(record.filled)))
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   identifier: ", hasidentifier(record) ? identifier(record) : "<missing>")
        println(io, "  description: ", hasdescription(record) ? description(record) : "<missing>")
          print(io, "     sequence: ", hassequence(record) ? truncate(sequence(String, record), 40) : "<missing>")
    else
        print(io, " <not filled>")
    end
end

function truncate(s::String, len::Integer)
    if length(s) > len
        return "$(String(collect(take(s, len - 1))))…"
    else
        return s
    end
end

function initialize!(record::Record)
    record.filled = 1:0
    record.identifier = 1:0
    record.description = 1:0
    record.sequence = 1:0
    return record
end

function Bio.isfilled(record::Record)
    return !isempty(record.filled)
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), p1, p2, n)
end


# Accessor functions
# ------------------

"""
    identifier(record::Record)::String

Get the sequence identifier of `record`.
"""
function identifier(record::Record)::String
    checkfilled(record)
    if !hasidentifier(record)
        missingerror(:identifier)
    end
    return String(record.data[record.identifier])
end

function hasidentifier(record)
    return isfilled(record) && !isempty(record.identifier)
end

function Bio.seqname(record::Record)
    return identifier(record)
end

function Bio.hasseqname(record::Record)
    return hasidentifier(record)
end

"""
    description(record::Record)::String

Get the description of `record`.
"""
function description(record::Record)::String
    checkfilled(record)
    if !hasdescription(record)
        missingerror(:description)
    end
    return String(record.data[record.description])
end

function hasdescription(record)
    return isfilled(record) && !isempty(record.description)
end

"""
    sequence(::Type{S}, record::Record, [part::UnitRange{Int}])::S

Get the sequence of `record`.

`S` can be either a subtype of `Bio.Seq.Sequence` or `String`.
If `part` argument is given, it returns the specified part of the sequence.
"""
function sequence{S<:Bio.Seq.Sequence}(::Type{S}, record::Record, part::UnitRange{Int}=1:endof(record.sequence))::S
    checkfilled(record)
    if !hassequence(record)
        missingerror(:sequence)
    end
    seqpart = record.sequence[part]
    return S(record.data, first(seqpart), last(seqpart))
end

function sequence(::Type{String}, record::Record, part::UnitRange{Int}=1:endof(record.sequence))::String
    checkfilled(record)
    if !hassequence(record)
        missingerror(:sequence)
    end
    return String(record.data[record.sequence[part]])
end

"""
    sequence(record::Record, [part::UnitRange{Int}])

Get the sequence of `record`.

This function infers the sequence type from the data. When it is wrong or
unreliable, use `sequence(::Type{S}, record::Record)`.  If `part` argument is
given, it returns the specified part of the sequence.
"""
function sequence(record::Record, part::UnitRange{Int}=1:endof(record.sequence))
    checkfilled(record)
    S = predict_seqtype(record.data, record.sequence)
    return sequence(S, record, part)
end

function hassequence(record::Record)
    # zero-length sequence may exist
    return isfilled(record)
end

function Bio.sequence(record::Record)
    return sequence(record)
end

function Bio.hassequence(record::Record)
    return hassequence(record)
end

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled FASTA record"))
    end
end

# Predict sequence type based on character frequencies in `seq[start:stop]`.
function predict_seqtype(seq::Vector{UInt8}, range)
    # count characters
    a = c = g = t = u = n = alpha = 0
    for i in range
        @inbounds x = seq[i]
        if x == 0x41 || x == 0x61
            a += 1
        elseif x == 0x43 || x == 0x63
            c += 1
        elseif x == 0x47 || x == 0x67
            g += 1
        elseif x == 0x54 || x == 0x74
            t += 1
        elseif x == 0x55 || x == 0x75
            u += 1
        elseif x == 0x4e || x == 0x6e
            n += 1
        end
        if 0x41 ≤ x ≤ 0x5a || 0x61 ≤ x ≤ 0x7a
            alpha += 1
            if alpha ≥ 300 && t + u > 0 && a + c + g + t + u + n == alpha
                # pretty sure that the sequence is either DNA or RNA
                break
            end
        end
    end

    # the threshold (= 0.95) is somewhat arbitrary
    if (a + c + g + t + u + n) / alpha > 0.95
        if t ≥ u
            return Bio.Seq.DNASequence
        else
            return Bio.Seq.RNASequence
        end
    else
        return Bio.Seq.AminoAcidSequence
    end
end
